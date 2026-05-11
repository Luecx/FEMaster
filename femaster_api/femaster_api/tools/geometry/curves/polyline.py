"""Polyline curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import clamp01, segment_lengths
from femaster_api.tools.geometry.vector import Point3, lerp, point3


@dataclass(frozen=True, slots=True)
class Polyline(Curve):
    points: tuple[Point3, ...] = ()

    def __init__(self, points: Iterable[Iterable[float]], *, seed: int = 1, tags: Iterable[str] = ()) -> None:
        items = tuple(point3(point) for point in points)
        if len(items) < 2:
            raise ValueError("polyline requires at least two points")
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "seed", int(seed))
        object.__setattr__(self, "points", items)
        self.__post_init__()

    def point_at(self, t: float) -> Point3:
        t = clamp01(t)
        lengths = segment_lengths(self.points)
        total = sum(lengths)
        if total <= 0.0:
            return self.points[0]
        target = t * total
        acc = 0.0
        for length, a, b in zip(lengths, self.points, self.points[1:]):
            if target <= acc + length:
                local = 0.0 if length <= 0.0 else (target - acc) / length
                return lerp(a, b, local)
            acc += length
        return self.points[-1]
