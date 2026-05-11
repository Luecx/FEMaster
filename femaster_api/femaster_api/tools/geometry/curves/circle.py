"""Circle curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, pi, sin
from typing import Iterable

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import clamp01
from femaster_api.tools.geometry.vector import Point3, add, distance, dot, plane_basis, point3, scale, sub


@dataclass(frozen=True, slots=True)
class Circle(Curve):
    center: Point3 = (0.0, 0.0, 0.0)
    radius: float = 1.0
    normal: Point3 = (0.0, 0.0, 1.0)
    reference: Point3 = (1.0, 0.0, 0.0)

    def __init__(
        self,
        center: Iterable[float],
        radius: float,
        *,
        normal: Iterable[float] = (0.0, 0.0, 1.0),
        reference: Iterable[float] | None = None,
        seed: int = 32,
        tags: Iterable[str] = (),
    ) -> None:
        if radius <= 0.0:
            raise ValueError("circle radius must be positive")
        u, _, n = plane_basis(normal, reference)
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "seed", int(seed))
        object.__setattr__(self, "center", point3(center))
        object.__setattr__(self, "radius", float(radius))
        object.__setattr__(self, "normal", n)
        object.__setattr__(self, "reference", u)
        self.__post_init__()

    @property
    def start(self) -> Point3:
        return self.point_at(0.0)

    @property
    def end(self) -> Point3:
        return self.point_at(1.0)

    def point_at(self, t: float) -> Point3:
        angle = 2.0 * pi * clamp01(t)
        u, v, _ = plane_basis(self.normal, self.reference)
        return add(self.center, add(scale(u, self.radius * cos(angle)), scale(v, self.radius * sin(angle))))

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        p = point3(point)
        _, _, n = plane_basis(self.normal, self.reference)
        rel = sub(p, self.center)
        return abs(dot(rel, n)) <= tolerance and abs(distance(p, self.center) - self.radius) <= tolerance
