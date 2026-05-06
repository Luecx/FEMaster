"""Orthotropic elasticity data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class OrthotropicElasticity:
    e1: float
    e2: float
    e3: float
    g23: float
    g13: float
    g12: float
    nu23: float
    nu13: float
    nu12: float
