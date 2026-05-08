"""Pressure load model object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.sets import SurfaceSet
from femaster_api.model.surfaces import SurfaceDefinition

SurfaceTarget = SurfaceDefinition | SurfaceSet


@dataclass(frozen=True, slots=True)
class PressureLoad:
    """Pressure applied normal to a surface or surface set."""

    target: SurfaceTarget
    pressure: float
    amplitude: Amplitude | None = None
