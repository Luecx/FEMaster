"""Amplitude data object."""

from __future__ import annotations

from dataclasses import dataclass

from .amplitude_interpolation import AmplitudeInterpolation


@dataclass(frozen=True, slots=True)
class Amplitude:
    name: str
    points: tuple[tuple[float, float], ...]
    interpolation: AmplitudeInterpolation = AmplitudeInterpolation.LINEAR
