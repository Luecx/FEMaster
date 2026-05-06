"""FEMaster serialization for generic fields."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.fields import FieldRepository


def write_fields(fields: FieldRepository) -> str:
    blocks: list[str] = []
    for field in fields.all():
        lines = [keyword("FIELD", NAME=field.name, TYPE=field.domain.value, COLS=field.cols, FILL=field.fill)]
        for row_id in sorted(field.values):
            lines.append(csv((row_id, *field.values[row_id])))
        blocks.append(block(lines))
    return "\n\n".join(block for block in blocks if block)
