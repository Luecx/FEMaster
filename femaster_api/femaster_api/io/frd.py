"""Reader for the nodal subset of CalculiX/CGX FRD written by FEMaster.

FEMaster's FRD writer emits NODE-domain result fields. This reader intentionally
focuses on that supported writer subset and maps it into the same Result/Field
objects used by the native RES reader.
"""

from __future__ import annotations

from pathlib import Path

from ..fields import Field, FieldDomain, FieldType
from ..results import Frame, LoadCase, Result


class FrdReader:
    """Parse FEMaster-generated ASCII FRD nodal result blocks."""

    def read(self, path: str | Path) -> Result:
        return self.parse(Path(path).read_text(encoding="utf-8", errors="replace"))

    def parse(self, text: str) -> Result:
        result = Result()
        loadcase = result.add(LoadCase(1))

        current_frame: Frame | None = None
        current_field: Field | None = None
        components: list[str] = []
        next_frame_id = 1

        for raw in text.splitlines():
            stripped = raw.strip()
            if not stripped:
                continue

            # FEMaster writes one 100CL control record for every output frame.
            # The exact fixed-width metadata is intentionally not exposed by the
            # public API; ordering is sufficient to reconstruct frame identity.
            if stripped.startswith("100CL"):
                current_frame = Frame(next_frame_id)
                next_frame_id += 1
                loadcase.add(current_frame)
                current_field = None
                components = []
                continue

            if stripped.startswith("-4"):
                if current_frame is None:
                    current_frame = loadcase.add(Frame(next_frame_id))
                    next_frame_id += 1

                tokens = stripped.split()
                name = tokens[1] if len(tokens) > 1 else "FIELD"
                current_field = Field(
                    name,
                    FieldDomain.NODE,
                    (),
                    type=FieldType.from_name(name),
                )
                current_frame.add(current_field)
                components = []
                continue

            if stripped.startswith("-5") and current_field is not None:
                tokens = stripped.split()
                if len(tokens) > 1:
                    components.append(tokens[1])
                    current_field.components = tuple(components)
                continue

            if stripped.startswith("-1") and current_field is not None:
                tokens = stripped.split()
                if len(tokens) < 3:
                    continue
                node_id = int(tokens[1])
                values = [float(value.replace("D", "E")) for value in tokens[2:]]

                # FRD -5 records may contain derived display components such as
                # vector norms which are not present in the -1 data record. Keep
                # only component names that correspond to actual values.
                if current_field.components and len(current_field.components) != len(values):
                    current_field.components = current_field.components[:len(values)]
                current_field.set(node_id, values)
                continue

            if stripped.startswith("-3"):
                current_field = None
                components = []

        return result
