"""Interpolation modes for amplitude objects."""

from __future__ import annotations

from enum import Enum


class AmplitudeInterpolation(Enum):
    LINEAR = "LINEAR"
    STEP = "STEP"
    NEAREST = "NEAREST"
