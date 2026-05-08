"""FEMaster serialization for generic fields."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.elements import Element
from femaster_api.model.fields import FieldDomain, FieldRepository
from femaster_api.model.nodes import Node


def write_fields(fields: FieldRepository, context: ExportContext) -> str:
    blocks: list[str] = []
    for field in fields:
        lines = [keyword("FIELD", NAME=field.name, TYPE=field.domain.value, COLS=field.cols, FILL=field.fill)]
        rows = sorted(field.values, key=lambda key: _field_sort_key(field.domain, key, context))
        for row_key in rows:
            lines.append(csv((*_field_key(row_key, context), *field.values[row_key])))
        blocks.append(block(lines))
    return "\n\n".join(block for block in blocks if block)


def _field_key(key: object, context: ExportContext) -> tuple[object, ...]:
    if isinstance(key, Node):
        return (context.node_id(key),)
    if isinstance(key, Element):
        return (context.element_id(key),)
    if isinstance(key, tuple):
        values: list[object] = []
        for item in key:
            if isinstance(item, Node):
                values.append(context.node_id(item))
            elif isinstance(item, Element):
                values.append(context.element_id(item))
            else:
                values.append(item)
        return tuple(values)
    return (key,)


def _field_sort_key(domain: FieldDomain, key: object, context: ExportContext) -> tuple[object, ...]:
    values = _field_key(key, context)
    return (domain.value, *((type(value).__name__, value) for value in values))
