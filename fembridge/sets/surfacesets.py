from __future__ import annotations
from typing import List, Optional, Iterable
from .surfaceset import SurfaceSet


class SurfaceSets:
    """
    Container for multiple SurfaceSet objects.
    Responsible for concatenating their FEMaster blocks.
    """

    def __init__(self) -> None:
        self._items: List[Optional[SurfaceSet]] = []

    def add(self, s: SurfaceSet) -> SurfaceSet:
        self._items.append(s)
        return s

    def extend(self, items: Iterable[SurfaceSet]) -> None:
        for s in items:
            self._items.append(s)

    def __getitem__(self, idx: int) -> SurfaceSet:
        val = self._items[idx]
        if val is None:
            raise KeyError(f"SurfaceSet index {idx} is empty.")
        return val

    def __len__(self) -> int:
        return len(self._items)

    # ------------- serialization -------------

    def to_femaster(self) -> str:
        blocks: List[str] = []
        for s in self._items:
            if s is None:
                continue
            blocks.append(s.to_femaster())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("SurfaceSets.to_asami not implemented yet.")
