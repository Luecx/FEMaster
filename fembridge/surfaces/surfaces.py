from __future__ import annotations
from typing import List, Iterable

from .surface import Surface

class Surfaces:
    """
    Kleiner interner Container für einzelne Surface-Objekte.
    Druck/Serialisierung erfolgt ausschließlich über SurfaceSet.
    """

    def __init__(self) -> None:
        self._items: List[Surface] = []

    def add(self, surface: Surface) -> Surface:
        self._items.append(surface)
        return surface

    def extend(self, items: Iterable[Surface]) -> None:
        self._items.extend(items)

    def __getitem__(self, idx: int) -> Surface:
        return self._items[idx]

    def __len__(self) -> int:
        return len(self._items)

    def items(self) -> List[Surface]:
        return list(self._items)

    def to_femaster(self) -> str:
        return ""

    def to_asami(self) -> str:
        raise NotImplementedError("Surface.to_asami not implemented yet.")
