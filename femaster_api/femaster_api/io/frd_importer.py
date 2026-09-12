"""Importer for FEMaster-generated ASCII FRD nodal results."""

from __future__ import annotations

from pathlib import Path

from ..model import Field, FieldDomain, FieldType, Frame, LoadCase, Result


class FrdImporter:
    """Import the nodal FRD subset written by FEMaster."""

    def import_file(self, path: str | Path) -> Result:
        return self.import_text(
            Path(path).read_text(encoding="utf-8", errors="replace")
        )

    def import_text(self, text: str) -> Result:
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
                values = [
                    float(value.replace("D", "E"))
                    for value in tokens[2:]
                ]

                if (
                    current_field.components
                    and len(current_field.components) != len(values)
                ):
                    current_field.components = current_field.components[:len(values)]

                current_field.set(node_id, values)
                continue

            if stripped.startswith("-3"):
                current_field = None
                components = []

        return result
