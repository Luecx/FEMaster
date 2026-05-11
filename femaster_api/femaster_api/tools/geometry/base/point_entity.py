"""0D geometry entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.vector import Point3, distance, point3

from .geometry_entity import GeometryEntity


@dataclass(frozen=True, slots=True)
class PointEntity(GeometryEntity):
    point: Point3 = (0.0, 0.0, 0.0)

    def __init__(self, point: Iterable[float], tags: Iterable[str] = ()) -> None:
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "point", point3(point))

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        return distance(self.point, point3(point)) <= tolerance
