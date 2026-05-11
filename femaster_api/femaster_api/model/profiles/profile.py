"""Beam profile data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Profile:
    name: str
    area: float
    iy: float
    iz: float
    j: float
    iyz: float = 0.0
    ey: float = 0.0
    ez: float = 0.0
    refy: float = 0.0
    refz: float = 0.0
