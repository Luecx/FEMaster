"""Pressure load data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads.amplitude import Amplitude

Target = str | int | object


@dataclass(frozen=True, slots=True)
class PressureLoad:
    target: Target
    pressure: float
    amplitude: Amplitude | None = None
