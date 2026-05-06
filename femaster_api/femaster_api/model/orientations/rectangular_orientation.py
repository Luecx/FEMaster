"""Rectangular orientation data object."""

from __future__ import annotations

from dataclasses import dataclass

Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class RectangularOrientation:
    name: str
    x_axis: Vector3
    y_axis: Vector3 | None = None
    z_axis: Vector3 | None = None
