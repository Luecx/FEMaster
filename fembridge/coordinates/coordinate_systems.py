from __future__ import annotations

from typing import List

from .coordinate_system_base import CoordinateSystem


class CoordinateSystems:
    """Container for coordinate system definitions."""

    def __init__(self) -> None:
        self._items: List[CoordinateSystem] = []

    def add(self, system: CoordinateSystem) -> CoordinateSystem:
        self._items.append(system)
        return system

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, idx: int) -> CoordinateSystem:
        return self._items[idx]

    def to_femaster(self) -> str:
        blocks: list[str] = []
        for system in self._items:
            blocks.append(system.to_femaster())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("CoordinateSystems.to_asami not implemented yet.")
