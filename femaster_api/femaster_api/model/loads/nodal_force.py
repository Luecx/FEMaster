"""Nodal force data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.orientations import Orientation

Target = str | int | object
Vector6 = tuple[float, float, float, float, float, float]


@dataclass(frozen=True, slots=True)
class NodalForce:
    target: Target
    values: Vector6
    orientation: Orientation | None = None
    amplitude: Amplitude | None = None
