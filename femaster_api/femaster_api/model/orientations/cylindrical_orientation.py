"""Cylindrical orientation data object."""

from __future__ import annotations

from dataclasses import dataclass

Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class CylindricalOrientation:
    name: str
    origin: Vector3
    axis: Vector3
    reference: Vector3
