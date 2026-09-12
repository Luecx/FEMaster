"""Interpolation modes for named FEMaster amplitudes.

The enum mirrors the native amplitude interpolation tokens and is kept in a
separate module so the ``Amplitude`` data object remains focused on sample
storage and export.
"""

from __future__ import annotations

from enum import Enum


class AmplitudeInterpolation(Enum):
    """Interpolation mode applied between successive amplitude samples."""

    LINEAR = "LINEAR"
    STEP = "STEP"
