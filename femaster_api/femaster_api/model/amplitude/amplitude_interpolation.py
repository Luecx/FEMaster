"""Amplitude interpolation modes."""

from enum import Enum


class AmplitudeInterpolation(Enum):
    """Interpolation mode between amplitude samples."""

    LINEAR = "LINEAR"
    STEP = "STEP"
