from __future__ import annotations
from typing import List
from .profile_base import Profile


class Profiles:
    """
    Container for Profile objects.

    - No IDs; profiles are referenced by their .name.
    - to_femaster() emits each profile block once (de-duped by name),
      in the order of first appearance.
    """
    def __init__(self) -> None:
        self._items: List[Profile] = []

    def add(self, profile: Profile) -> Profile:
        self._items.append(profile)
        return profile

    def __getitem__(self, idx: int) -> Profile:
        return self._items[idx]

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        blocks: List[str] = []
        seen: set[str] = set()
        for p in self._items:
            if p.name in seen:
                continue
            seen.add(p.name)
            blocks.append(p.to_femaster().strip())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("Profiles.to_asami not implemented yet.")
