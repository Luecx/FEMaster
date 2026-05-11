"""Joined curve entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.tools.geometry.base import Curve
from femaster_api.tools.geometry.base.curve_math import clamp01, curve_length
from femaster_api.tools.geometry.vector import Point3, distance


@dataclass(frozen=True, slots=True)
class JoinedCurve(Curve):
    curves: tuple[Curve, ...] = ()

    def __init__(self, curves: Iterable[Curve], *, tags: Iterable[str] = ()) -> None:
        items = tuple(curves)
        if len(items) < 1:
            raise ValueError("joined curve requires at least one curve")
        for left, right in zip(items, items[1:]):
            if distance(left.end, right.start) > 1.0e-6:
                raise ValueError("curves do not connect head-to-tail")
        seed = sum(curve.subdivisions for curve in items)
        inherited_tags = tuple(tag for curve in items for tag in curve.tags)
        object.__setattr__(self, "tags", tuple(dict.fromkeys((*inherited_tags, *tuple(tags)))))
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "curves", items)
        self.__post_init__()

    def point_at(self, t: float) -> Point3:
        t = clamp01(t)
        lengths = tuple(curve_length(curve) for curve in self.curves)
        total = sum(lengths)
        if total <= 0.0:
            return self.curves[0].start
        target = t * total
        acc = 0.0
        for length, curve in zip(lengths, self.curves):
            if target <= acc + length:
                local = 0.0 if length <= 0.0 else (target - acc) / length
                return curve.point_at(local)
            acc += length
        return self.curves[-1].end
