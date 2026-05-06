"""Material data object."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .abd_elasticity import ABDElasticity
from .generalized_isotropic_elasticity import GeneralizedIsotropicElasticity
from .isotropic_elasticity import IsotropicElasticity
from .orthotropic_elasticity import OrthotropicElasticity

Elasticity = IsotropicElasticity | GeneralizedIsotropicElasticity | OrthotropicElasticity | ABDElasticity


@dataclass(frozen=True, slots=True)
class Material:
    name: str
    elasticity: Elasticity | None = None
    density: float | None = None
    thermal_expansion: float | None = None

    def set_elasticity(self, elasticity: Elasticity) -> "Material":
        if not isinstance(
            elasticity,
            (IsotropicElasticity, GeneralizedIsotropicElasticity, OrthotropicElasticity, ABDElasticity),
        ):
            raise TypeError("elasticity must be an elasticity object")
        return replace(self, elasticity=elasticity)

    def set_isotropic_elasticity(
        self,
        youngs_modulus: float | None = None,
        poisson_ratio: float | None = None,
        *,
        E: float | None = None,
        nu: float | None = None,
    ) -> "Material":
        e_value = youngs_modulus if youngs_modulus is not None else E
        nu_value = poisson_ratio if poisson_ratio is not None else nu
        if e_value is None or nu_value is None:
            raise ValueError("isotropic elasticity requires youngs_modulus/E and poisson_ratio/nu")
        return self.set_elasticity(IsotropicElasticity(float(e_value), float(nu_value)))

    def set_generalized_isotropic_elasticity(
        self,
        youngs_modulus: float,
        poisson_ratio: float,
        shear_modulus: float,
    ) -> "Material":
        return self.set_elasticity(
            GeneralizedIsotropicElasticity(
                float(youngs_modulus),
                float(poisson_ratio),
                float(shear_modulus),
            )
        )

    def set_orthotropic_elasticity(
        self,
        *,
        e1: float,
        e2: float,
        e3: float,
        g23: float,
        g13: float,
        g12: float,
        nu23: float,
        nu13: float,
        nu12: float,
    ) -> "Material":
        return self.set_elasticity(
            OrthotropicElasticity(
                float(e1),
                float(e2),
                float(e3),
                float(g23),
                float(g13),
                float(g12),
                float(nu23),
                float(nu13),
                float(nu12),
            )
        )

    def set_abd_elasticity(self, values: tuple[float, ...] | list[float]) -> "Material":
        return self.set_elasticity(ABDElasticity(tuple(float(value) for value in values)))

    def set_density(self, density: float | None) -> "Material":
        return replace(self, density=None if density is None else float(density))

    def set_thermal_expansion(self, thermal_expansion: float | None) -> "Material":
        return replace(self, thermal_expansion=None if thermal_expansion is None else float(thermal_expansion))
