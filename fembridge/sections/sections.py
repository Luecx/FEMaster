from __future__ import annotations
from typing import List
from .section_base import Section


class Sections:
    """Container for all Section objects (beam, shell, solid, etc.)."""
    def __init__(self) -> None:
        self._items: List[Section] = []

    def add(self, section: Section) -> Section:
        """Add a Section (no IDs or names needed)."""
        self._items.append(section)
        return section

    def __getitem__(self, idx: int) -> Section:
        return self._items[idx]

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        """Join all section definitions into one FEMaster block string."""
        blocks: List[str] = []
        for s in self._items:
            block = s.to_femaster().strip()
            if block:
                blocks.append(block)
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("Sections.to_asami not implemented yet.")
