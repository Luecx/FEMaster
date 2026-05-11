"""Base 1D geometry entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base.curve_math import point_on_segment
from femaster_api.tools.geometry.vector import Point3, point3

from .geometry_entity import GeometryEntity


@dataclass(frozen=True, slots=True)
class Curve(GeometryEntity):
    seed: int = 1

    def __post_init__(self) -> None:
        if self.seed < 1:
            raise ValueError("curve seed must be at least one")

    @property
    def subdivisions(self) -> int:
        return self.seed

    @property
    def start(self) -> Point3:
        return self.point_at(0.0)

    @property
    def end(self) -> Point3:
        return self.point_at(1.0)

    def point_at(self, t: float) -> Point3:
        raise NotImplementedError

    def sample(self, count: int | None = None) -> tuple[Point3, ...]:
        n = self.seed if count is None else max(1, int(count) - 1)
        return tuple(self.point_at(i / n) for i in range(n + 1))

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        target = point3(point)
        points = self.sample(max(2, self.seed * 8 + 1))
        return any(point_on_segment(target, a, b, tolerance) for a, b in zip(points, points[1:]))
