"""Generalized isotropic elasticity data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class GeneralizedIsotropicElasticity:
    youngs_modulus: float
    poisson_ratio: float
    shear_modulus: float
