"""ID extraction for object-backed sets."""

from __future__ import annotations

def unique_members(members: object) -> tuple[object, ...]:
    by_id: dict[int, object] = {}
    for member in members:
        by_id[id(member)] = member
    return tuple(by_id[key] for key in sorted(by_id))
