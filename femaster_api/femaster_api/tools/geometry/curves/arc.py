"""Circular arc curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from math import atan2, cos, pi, sin
from typing import Iterable

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import angle_between, clamp01
from femaster_api.tools.geometry.vector import Point3, add, distance, dot, plane_basis, point3, scale, sub


@dataclass(frozen=True, slots=True)
class Arc(Curve):
    center: Point3 = (0.0, 0.0, 0.0)
    radius: float = 1.0
    start_angle: float = 0.0
    end_angle: float = pi / 2.0
    normal: Point3 = (0.0, 0.0, 1.0)
    reference: Point3 = (1.0, 0.0, 0.0)

    def __init__(
        self,
        center: Iterable[float],
        radius: float,
        start_angle: float,
        end_angle: float,
        *,
        normal: Iterable[float] = (0.0, 0.0, 1.0),
        reference: Iterable[float] | None = None,
        seed: int = 8,
        tags: Iterable[str] = (),
    ) -> None:
        if radius <= 0.0:
            raise ValueError("arc radius must be positive")
        u, _, n = plane_basis(normal, reference)
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "seed", int(seed))
        object.__setattr__(self, "center", point3(center))
        object.__setattr__(self, "radius", float(radius))
        object.__setattr__(self, "start_angle", float(start_angle))
        object.__setattr__(self, "end_angle", float(end_angle))
        object.__setattr__(self, "normal", n)
        object.__setattr__(self, "reference", u)
        self.__post_init__()

    def point_at(self, t: float) -> Point3:
        angle = self.start_angle + clamp01(t) * (self.end_angle - self.start_angle)
        u, v, _ = plane_basis(self.normal, self.reference)
        return add(self.center, add(scale(u, self.radius * cos(angle)), scale(v, self.radius * sin(angle))))

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        p = point3(point)
        u, v, n = plane_basis(self.normal, self.reference)
        rel = sub(p, self.center)
        if abs(dot(rel, n)) > tolerance:
            return False
        radial = distance(p, self.center)
        if abs(radial - self.radius) > tolerance:
            return False
        angle = atan2(dot(rel, v), dot(rel, u))
        return angle_between(angle, self.start_angle, self.end_angle, tolerance)
