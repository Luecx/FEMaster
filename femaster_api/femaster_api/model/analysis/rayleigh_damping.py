"""Rayleigh damping data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class RayleighDamping:
    alpha: float
    beta: float
