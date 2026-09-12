"""Complete FEMaster result hierarchy and native result-file readers.

``Result`` is the public entry point for post-processing.  It owns one or more
``Solution`` objects; every solution owns ``Frame`` objects; every frame owns
``Field`` objects.  A frame has both a discrete ``id`` and an optional physical
``value``.  FRD supplies that value directly through its ``100CL`` record.  The
current native RES writer does not persist ``frame_value``; in that case the
model intentionally leaves ``Frame.value`` as ``None`` instead of inventing a
number.

Result reading lives here rather than in a separate ``io`` package.  The public
API is therefore format-oriented and discoverable from the result object itself:

    Result.read("job.res")
    Result.read_res("job.res")
    Result.read_frd("job.frd")

The parsers preserve semantic node/element identifiers, including qualified
instance-local identifiers such as ``"bolt.17"``.  Format-specific syntax is
normalized into the central ``FieldDomain`` and ``FieldType`` definitions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

from ..field.field import Field
from ..field.field_domain import FieldDomain
from ..field.field_type import FieldType
from ..field.typing import FieldIndex, FieldKey
from .result_frame import Frame
from .result_solution import Solution


class Result:
    """Complete format-independent result of one FEMaster execution."""

    def __init__(self) -> None:
        self.solutions: list[Solution] = []

    # ------------------------------------------------------------------
    # Public hierarchy access
    # ------------------------------------------------------------------

    def add(self, solution: Solution) -> Solution:
        """Append one solution while enforcing a unique numeric ID."""

        if any(item.id == solution.id for item in self.solutions):
            raise ValueError(f"duplicate solution id: {solution.id}")
        self.solutions.append(solution)
        return solution

    def solution(self, key: int | str = 1) -> Solution:
        """Return a solution by numeric ID or optional semantic name."""

        for solution in self.solutions:
            if isinstance(key, int) and solution.id == key:
                return solution
            if isinstance(key, str) and solution.name == key:
                return solution
        raise KeyError(f"unknown solution: {key}")

    def field(
        self,
        key: str | FieldType,
        *,
        solution: int | str = 1,
        frame: int | str = 0,
    ) -> Field:
        """Return one field directly from a selected solution and frame."""

        return self.solution(solution).field(key, frame=frame)

    # ------------------------------------------------------------------
    # Public file readers
    # ------------------------------------------------------------------

    @classmethod
    def read(cls, path: str | Path) -> "Result":
        """Read a supported result format selected by filename extension."""

        path = Path(path)
        suffix = path.suffix.lower()

        if suffix == ".res":
            return cls.read_res(path)
        if suffix == ".frd":
            return cls.read_frd(path)

        raise ValueError(
            f"unsupported result format: {path.suffix or '<none>'}"
        )

    @classmethod
    def read_res(cls, path: str | Path) -> "Result":
        """Read a native FEMaster text ``.res`` file."""

        text = Path(path).read_text(encoding="utf-8")
        return cls.read_res_text(text)

    @classmethod
    def read_res_text(cls, text: str) -> "Result":
        """Parse native FEMaster RES text into solutions, frames and fields."""

        result = cls()
        current_solution: Solution | None = None
        current_frame: Frame | None = None
        lines = iter(enumerate(text.splitlines(), start=1))

        for line_number, raw in lines:
            line = raw.strip()
            if not line or cls._is_comment(line):
                continue

            upper = line.upper()

            # ``LC`` starts one solver solution.  Current ResWriter emits only
            # the numeric id, but an optional trailing name is accepted so the
            # reader remains compatible with richer text output.
            if upper.startswith(("LC ", "LOADCASE ")):
                parts = line.split(maxsplit=2)
                solution_id = int(parts[1])
                solution_name = parts[2] if len(parts) > 2 else None
                current_solution = result.add(
                    Solution(solution_id, name=solution_name)
                )
                current_frame = None
                continue

            # FRAME is supported as an explicit forward-compatible record even
            # though the current native writer does not emit it.  Both
            # ``FRAME 2 0.5`` and ``FRAME, ID=2, VALUE=0.5`` are accepted.
            if upper.startswith("FRAME"):
                if current_solution is None:
                    current_solution = result.add(Solution(1))
                current_frame = cls._parse_res_frame(line)
                current_solution.add(current_frame)
                continue

            if upper.startswith("FIELD"):
                if current_solution is None:
                    current_solution = result.add(Solution(1))
                if current_frame is None:
                    current_frame = current_solution.add(Frame(0))

                field = cls._read_res_field(line, lines, line_number)
                current_frame.add(field)

        return result

    @classmethod
    def read_frd(cls, path: str | Path) -> "Result":
        """Read a FEMaster-generated ASCII CalculiX/CGX ``.frd`` file."""

        text = Path(path).read_text(
            encoding="utf-8",
            errors="replace",
        )
        return cls.read_frd_text(text)

    @classmethod
    def read_frd_text(cls, text: str) -> "Result":
        """Parse the nodal FRD result subset written by FEMaster."""

        result = cls()

        # Metadata records preceding every field block identify the solver step,
        # frame number and physical frame value.  They are retained until the
        # following -4 field-definition record creates the actual Field object.
        pending_step = 1
        pending_frame = 1
        pending_value: float | None = None
        pending_type: str | None = None

        current_frame: Frame | None = None
        current_field: Field | None = None
        current_node: FieldKey | None = None
        component_names: list[str] = []

        for raw in text.splitlines():
            stripped = raw.strip()
            if not stripped:
                continue

            # 1PSTEP contains block id, frame id and analysis-step id.
            if stripped.startswith("1PSTEP"):
                tokens = stripped.split()
                if len(tokens) >= 4:
                    pending_frame = int(tokens[-2])
                    pending_step = int(tokens[-1])
                continue

            # 100CL carries the physical frame value and increment type.  The
            # FEMaster writer uses split-friendly fixed-width fields:
            # 100CL, global id, value, row count, ictype, global frame, ...
            if stripped.startswith("100CL"):
                tokens = stripped.split()
                if len(tokens) >= 4:
                    pending_value = cls._float_token(tokens[2])

                if len(tokens) >= 5:
                    try:
                        ictype = int(tokens[4])
                    except ValueError:
                        ictype = 0
                    pending_type = {
                        0: "STATIC",
                        1: "DYNAMIC",
                        2: "MODAL",
                        4: "BUCKLING",
                    }.get(ictype)

                solution = cls._get_or_create_solution(
                    result,
                    pending_step,
                    type=pending_type,
                )
                current_frame = cls._get_or_create_frame(
                    solution,
                    pending_frame,
                    pending_value,
                )
                continue

            # -4 starts the field definition belonging to the preceding frame
            # metadata.  If metadata is absent, create a conservative default.
            if stripped.startswith("-4"):
                solution = cls._get_or_create_solution(
                    result,
                    pending_step,
                    type=pending_type,
                )
                if current_frame is None:
                    current_frame = cls._get_or_create_frame(
                        solution,
                        pending_frame,
                        pending_value,
                    )

                tokens = stripped.split()
                name = tokens[1] if len(tokens) > 1 else "FIELD"
                current_field = Field(
                    name,
                    FieldDomain.NODE,
                    (),
                    type=FieldType.from_name(name),
                )
                current_frame.add(current_field)
                current_node = None
                component_names = []
                continue

            # -5 declares one component, including derived display components.
            if stripped.startswith("-5") and current_field is not None:
                tokens = stripped.split()
                if len(tokens) > 1:
                    component_names.append(tokens[1])
                    current_field.components = tuple(component_names)
                continue

            # -1 starts one nodal result row.  Mesh -1 records are ignored
            # because no current result field exists while geometry is parsed.
            if stripped.startswith("-1") and current_field is not None:
                tokens = stripped.split()
                if len(tokens) < 3:
                    continue
                current_node = cls._semantic_id(tokens[1])
                values = [
                    cls._float_token(token)
                    for token in tokens[2:]
                ]
                current_field.set(current_node, values)
                continue

            # -2 continues the current node when a result has more than six
            # components.  Append rather than replacing the existing row.
            if stripped.startswith("-2") and current_field is not None:
                if current_node is None:
                    continue
                tokens = stripped.split()
                continuation = tuple(
                    cls._float_token(token)
                    for token in tokens[1:]
                )
                previous = current_field.values[current_node]
                current_field.values[current_node] = previous + continuation
                continue

            # -3 closes the numerical field block.
            if stripped.startswith("-3") and current_field is not None:
                current_field = None
                current_node = None
                component_names = []
                # A following field may share the same frame; the next 100CL
                # record will re-resolve it without creating a duplicate frame.
                continue

        return result

    # ------------------------------------------------------------------
    # Native RES parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_res_frame(header: str) -> Frame:
        normalized = header.replace(",", " ")
        tokens = normalized.split()

        keys: dict[str, str] = {}
        positional: list[str] = []
        for token in tokens[1:]:
            if "=" in token:
                key, value = token.split("=", 1)
                keys[key.upper()] = value
            else:
                positional.append(token)

        frame_id = int(keys.get("ID", positional[0] if positional else "0"))

        raw_value = keys.get("VALUE")
        if raw_value is None and len(positional) > 1:
            raw_value = positional[1]

        value = None if raw_value is None else float(raw_value)
        name = keys.get("NAME")
        return Frame(frame_id, value=value, name=name)

    @classmethod
    def _read_res_field(
        cls,
        header: str,
        lines: Iterator[tuple[int, str]],
        header_line: int,
    ) -> Field:
        command, keys = cls._parse_res_field_header(header)

        rows = int(keys.get("ROWS", "0"))
        name = keys.get("NAME", command)
        index_cols = int(keys.get("INDEX_COLS", "0"))

        if "VALUE_COLS" in keys:
            value_cols = int(keys["VALUE_COLS"])
        elif "COLS" in keys:
            value_cols = int(keys["COLS"])
        else:
            value_cols = 0

        components = tuple(
            item.strip()
            for item in keys.get("COMPONENTS", "").split(";")
            if item.strip()
        )
        if not components:
            components = tuple(
                f"C{index + 1}"
                for index in range(value_cols)
            )

        field = Field(
            name,
            cls._res_field_domain(keys, index_cols),
            components,
            type=FieldType.from_name(name),
        )

        read_rows = 0
        for line_number, raw in lines:
            line = raw.strip()
            if not line or cls._is_comment(line):
                continue
            if line.upper().startswith("END FIELD"):
                break

            tokens = line.replace(",", " ").split()
            expected = index_cols + value_cols
            if len(tokens) < expected:
                raise ValueError(
                    f"RES field {name!r} at line {header_line} expects "
                    f"{expected} columns, got {len(tokens)} at line {line_number}"
                )

            key = cls._res_field_key(tokens[:index_cols], read_rows)
            values = tuple(
                cls._float_token(token)
                for token in tokens[
                    index_cols:index_cols + value_cols
                ]
            )
            field.set(key, values)
            read_rows += 1

            if rows and read_rows >= rows:
                break

        if rows and read_rows != rows:
            raise ValueError(
                f"RES field {name!r} expected {rows} rows, got {read_rows}"
            )

        return field

    @staticmethod
    def _parse_res_field_header(
        header: str,
    ) -> tuple[str, dict[str, str]]:
        normalized = header.replace(",", " ")
        parts = normalized.split()

        if len(parts) >= 2 and "=" not in parts[1]:
            command = parts[1]
            tokens = parts[2:]
        else:
            command = "FIELD"
            tokens = parts[1:]

        keys: dict[str, str] = {}
        for token in tokens:
            if "=" in token:
                key, value = token.split("=", 1)
                keys[key.upper()] = value

        return command, keys

    @staticmethod
    def _res_field_domain(
        keys: dict[str, str],
        index_cols: int,
    ) -> FieldDomain:
        raw = keys.get("TYPE") or keys.get("DOMAIN")
        if raw:
            normalized = raw.replace("_", "").upper()
            if normalized == "IP":
                return FieldDomain.ELEMENT_IP

            for domain in FieldDomain:
                candidates = {
                    domain.value.replace("_", "").upper(),
                    domain.name.replace("_", "").upper(),
                }
                if normalized in candidates:
                    return domain

        if index_cols > 0:
            return FieldDomain.ELEMENT_NODAL
        return FieldDomain.UNKNOWN

    @classmethod
    def _res_field_key(
        cls,
        tokens: list[str],
        fallback: int,
    ) -> FieldKey:
        if not tokens:
            return fallback

        values = tuple(cls._semantic_id(token) for token in tokens)
        return values[0] if len(values) == 1 else values

    # ------------------------------------------------------------------
    # Shared result parsing helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _semantic_id(token: str) -> FieldIndex:
        """Keep bare integer IDs numeric and qualified instance IDs textual."""

        try:
            return int(token)
        except ValueError:
            return token

    @staticmethod
    def _float_token(token: str) -> float:
        """Parse ordinary or Fortran-style floating-point result text."""

        return float(token.replace("D", "E").replace("d", "e"))

    @staticmethod
    def _is_comment(line: str) -> bool:
        return line.startswith(("#", "!", "//", "**"))

    @staticmethod
    def _get_or_create_solution(
        result: "Result",
        id: int,
        *,
        type: str | None = None,
    ) -> Solution:
        for solution in result.solutions:
            if solution.id == id:
                if solution.type is None and type is not None:
                    solution.type = type
                return solution

        return result.add(Solution(id, type=type))

    @staticmethod
    def _get_or_create_frame(
        solution: Solution,
        id: int,
        value: float | None,
    ) -> Frame:
        for frame in solution.frames:
            if frame.id == id:
                if frame.value is None and value is not None:
                    frame.value = value
                return frame

        return solution.add(Frame(id, value=value))

    def __len__(self) -> int:
        return len(self.solutions)
