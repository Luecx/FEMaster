"""Inertial load data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import EntitySet

Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class InertialLoad:
    target: EntitySet
    center: Vector3
    center_acceleration: Vector3
    omega: Vector3
    alpha: Vector3
    consider_point_masses: bool = False
