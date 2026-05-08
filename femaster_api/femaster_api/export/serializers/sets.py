"""FEMaster serialization for object-backed entity sets."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.sets import (
    ElementSetRepository,
    EntitySet,
    EntityType,
    NodeSetRepository,
    SurfaceSetRepository,
)


def write_sets(
    node_sets: NodeSetRepository,
    element_sets: ElementSetRepository,
    surface_sets: SurfaceSetRepository,
    context: ExportContext,
) -> str:
    items = (*tuple(node_sets), *tuple(element_sets), *tuple(surface_sets))
    blocks = [_write_set(item, context) for item in items]
    return "\n\n".join(block for block in blocks if block)


def _write_set(item: EntitySet, context: ExportContext) -> str:
    ids = _ids(item, context)
    if not ids:
        return ""
    command = {
        EntityType.NODE: "NSET",
        EntityType.ELEMENT: "ELSET",
        EntityType.SURFACE: "SFSET",
    }[item.entity_type]
    key = command
    if _is_regular_range(ids):
        step = ids[1] - ids[0] if len(ids) > 1 else 1
        return block([keyword(command, **{key: item.name, "GENERATE": True}), csv((ids[0], ids[-1], step))])

    lines = [keyword(command, **{key: item.name})]
    for start in range(0, len(ids), 16):
        lines.append(csv(ids[start : start + 16]))
    return block(lines)


def _is_regular_range(ids: tuple[int, ...]) -> bool:
    if len(ids) <= 2:
        return True
    step = ids[1] - ids[0]
    return all((right - left) == step for left, right in zip(ids, ids[1:]))


def _ids(item: EntitySet, context: ExportContext) -> tuple[int, ...]:
    if item.entity_type == EntityType.NODE:
        return tuple(context.node_id(member) for member in item.members)
    if item.entity_type == EntityType.ELEMENT:
        return tuple(context.element_id(member) for member in item.members)
    if item.entity_type == EntityType.SURFACE:
        return tuple(context.surface_id(member) for member in item.members)
    raise TypeError(f"unsupported set type: {item.entity_type}")
