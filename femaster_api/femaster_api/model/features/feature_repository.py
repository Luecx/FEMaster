"""Repository for model features."""

from __future__ import annotations

from typing import Iterator

from .point_mass import PointMass


class FeatureRepository:
    """Repository for point masses and future model features."""

    def __init__(self) -> None:
        self._point_masses: list[PointMass] = []

    def add(self, point_mass: PointMass) -> PointMass:
        if not isinstance(point_mass, PointMass):
            raise TypeError("point_mass must be a PointMass")
        self._point_masses.append(point_mass)
        return point_mass

    def point_masses(self) -> tuple[PointMass, ...]:
        return tuple(self._point_masses)

    def __iter__(self) -> Iterator[PointMass]:
        return iter(self._point_masses)

    def __getitem__(self, index: int | slice) -> PointMass | tuple[PointMass, ...]:
        if isinstance(index, slice):
            return tuple(self._point_masses[index])
        return self._point_masses[index]

    def __len__(self) -> int:
        return len(self._point_masses)
