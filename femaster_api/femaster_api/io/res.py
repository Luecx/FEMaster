"""Reader for FEMaster text result (RES) files."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

from ..fields import Field, FieldDomain, FieldType
from ..results import Frame, LoadCase, Result


class ResReader:
    """Parse FEMaster RES text into the common result hierarchy."""

    def read(self, path: str | Path) -> Result:
        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> Result:
        result = Result()
        current_loadcase: LoadCase | None = None
        current_frame: Frame | None = None
        lines = iter(enumerate(text.splitlines(), start=1))

        for _, raw in lines:
            line = raw.strip()
            if not line or self._is_comment(line):
                continue

            upper = line.upper()

            if upper.startswith(("LC ", "LOADCASE ")):
                parts = line.split(maxsplit=2)
                current_loadcase = LoadCase(
                    int(parts[1]),
                    parts[2] if len(parts) > 2 else None,
                )
                current_frame = None
                result.add(current_loadcase)
                continue

            if upper.startswith("FRAME"):
                if current_loadcase is None:
                    current_loadcase = result.add(LoadCase(1))
                current_frame = self._parse_frame(line)
                current_loadcase.add(current_frame)
                continue

            if upper.startswith("FIELD"):
                if current_loadcase is None:
                    current_loadcase = result.add(LoadCase(1))
                if current_frame is None:
                    current_frame = current_loadcase.add(Frame())
                current_frame.add(self._read_field(line, lines))

        return result

    @staticmethod
    def _parse_frame(header: str) -> Frame:
        parts = header.split(maxsplit=2)
        if len(parts) == 1:
            return Frame()
        try:
            return Frame(int(parts[1]), parts[2] if len(parts) > 2 else None)
        except ValueError:
            return Frame(0, header.partition(" ")[2] or None)

    def _read_field(self, header: str, lines: Iterator[tuple[int, str]]) -> Field:
        command, keys = self._parse_field_header(header)
        rows       = int(keys.get("ROWS", "0"))
        name       = keys.get("NAME", command)
        index_cols = int(keys.get("INDEX_COLS", "0"))

        if "VALUE_COLS" in keys:
            value_cols = int(keys["VALUE_COLS"])
        elif "COLS" in keys:
            value_cols = int(keys["COLS"]) - index_cols
        else:
            value_cols = 0

        components = tuple(
            item.strip()
            for item in keys.get("COMPONENTS", "").split(";")
            if item.strip()
        )
        if not components:
            components = tuple(f"C{i + 1}" for i in range(value_cols))

        field = Field(
            name,
            self._field_domain(keys, index_cols),
            components,
            type=FieldType.from_name(name),
        )

        read_rows = 0
        for _, raw in lines:
            line = raw.strip()
            if not line or self._is_comment(line):
                continue
            if line.upper().startswith("END FIELD"):
                break

            values = tuple(float(token) for token in line.replace(",", " ").split())
            if value_cols == 0:
                value_cols = max(0, len(values) - index_cols)
            key = self._field_key(values, index_cols, read_rows)
            field.set(key, values[index_cols:index_cols + value_cols])
            read_rows += 1

            if rows and read_rows >= rows:
                break

        return field

    @staticmethod
    def _field_key(
        values: tuple[float, ...],
        index_cols: int,
        fallback: int,
    ) -> int | tuple[int, ...]:
        if index_cols == 0:
            return fallback
        index = tuple(int(value) for value in values[:index_cols])
        return index[0] if len(index) == 1 else index

    @staticmethod
    def _field_domain(keys: dict[str, str], index_cols: int) -> FieldDomain:
        raw = keys.get("TYPE") or keys.get("DOMAIN")
        if raw:
            normalized = raw.replace("_", "").upper()
            if normalized == "IP":
                return FieldDomain.ELEMENT_IP
            for domain in FieldDomain:
                if normalized in {
                    domain.value.replace("_", "").upper(),
                    domain.name.replace("_", "").upper(),
                }:
                    return domain
        if index_cols > 0:
            return FieldDomain.ELEMENT_NODAL
        return FieldDomain.UNKNOWN

    @staticmethod
    def _parse_field_header(header: str) -> tuple[str, dict[str, str]]:
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
    def _is_comment(line: str) -> bool:
        return line.startswith(("#", "!", "//", "**"))
