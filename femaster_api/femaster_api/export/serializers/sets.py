"""FEMaster serialization for entity sets."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.sets import EntitySet, EntityType, SetRepository


def write_sets(sets: SetRepository) -> str:
    blocks = [_write_set(item) for item in sets.all()]
    return "\n\n".join(block for block in blocks if block)


def _write_set(item: EntitySet) -> str:
    ids = _ids(item)
    if not ids:
        return ""
    command = {
        EntityType.NODE: "NSET",
        EntityType.ELEMENT: "ELSET",
        EntityType.SURFACE: "SFSET",
    }[item.entity_type]
    key = command
    if item.generated and _is_regular_range(ids):
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


def _ids(item: EntitySet) -> tuple[int, ...]:
    return tuple(member.id for member in item.members)
