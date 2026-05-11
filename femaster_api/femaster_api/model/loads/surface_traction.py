"""Surface traction data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.orientations import Orientation

Target = str | int | object
Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class SurfaceTraction:
    target: Target
    values: Vector3
    orientation: Orientation | None = None
    amplitude: Amplitude | None = None
