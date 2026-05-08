"""Small helpers for object-backed sets."""

from __future__ import annotations

def unique_members(members: object) -> tuple[object, ...]:
    """Return unique members by object identity while preserving input order."""

    seen: set[int] = set()
    unique: list[object] = []
    for member in members:
        marker = id(member)
        if marker not in seen:
            seen.add(marker)
            unique.append(member)
    return tuple(unique)
