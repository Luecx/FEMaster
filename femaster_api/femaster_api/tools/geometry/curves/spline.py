"""Spline curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import catmull_rom, clamp01
from femaster_api.tools.geometry.curves.polyline import Polyline
from femaster_api.tools.geometry.vector import Point3, add, point3, scale

Interpolation = Literal["linear", "bezier", "catmull_rom"]


@dataclass(frozen=True, slots=True)
class Spline(Curve):
    control_points: tuple[Point3, ...] = ()
    interpolation: Interpolation = "catmull_rom"

    def __init__(
        self,
        control_points: Iterable[Iterable[float]],
        *,
        interpolation: Interpolation = "catmull_rom",
        seed: int = 16,
        tags: Iterable[str] = (),
    ) -> None:
        points = tuple(point3(point) for point in control_points)
        if len(points) < 2:
            raise ValueError("spline requires at least two control points")
        if interpolation == "bezier" and len(points) != 4:
            raise ValueError("bezier spline currently requires exactly four control points")
        if interpolation not in ("linear", "bezier", "catmull_rom"):
            raise ValueError(f"unknown interpolation: {interpolation}")
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "seed", int(seed))
        object.__setattr__(self, "control_points", points)
        object.__setattr__(self, "interpolation", interpolation)
        self.__post_init__()

    def point_at(self, t: float) -> Point3:
        t = clamp01(t)
        if self.interpolation == "linear":
            return Polyline(self.control_points, seed=self.seed).point_at(t)
        if self.interpolation == "bezier":
            p0, p1, p2, p3 = self.control_points
            return add(
                add(scale(p0, (1 - t) ** 3), scale(p1, 3 * (1 - t) ** 2 * t)),
                add(scale(p2, 3 * (1 - t) * t**2), scale(p3, t**3)),
            )
        return catmull_rom(self.control_points, t)
