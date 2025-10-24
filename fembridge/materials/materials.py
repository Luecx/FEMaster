
from __future__ import annotations
from typing import List, Optional
from .material_base import Material

class Materials:
    """Simple container for Material objects (materials referenced by name)."""
    def __init__(self) -> None:
        self._items: List[Optional[Material]] = []

    def add(self, mat: Material) -> Material:
        self._items.append(mat)
        return mat

    def to_femaster(self) -> str:
        blocks = []
        for m in self._items:
            if m is None:
                continue
            blocks.append(m.to_femaster())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("Materials.to_asami not implemented yet.")
