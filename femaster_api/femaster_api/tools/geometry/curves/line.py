"""Line curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import clamp01
from femaster_api.tools.geometry.vector import Point3, lerp, point3


@dataclass(frozen=True, slots=True)
class Line(Curve):
    a: Point3 = (0.0, 0.0, 0.0)
    b: Point3 = (1.0, 0.0, 0.0)

    def __init__(self, a: Iterable[float], b: Iterable[float], *, seed: int = 1, tags: Iterable[str] = ()) -> None:
        object.__setattr__(self, "tags", tuple(tags))
        object.__setattr__(self, "seed", int(seed))
        object.__setattr__(self, "a", point3(a))
        object.__setattr__(self, "b", point3(b))
        self.__post_init__()

    def point_at(self, t: float) -> Point3:
        return lerp(self.a, self.b, clamp01(t))
