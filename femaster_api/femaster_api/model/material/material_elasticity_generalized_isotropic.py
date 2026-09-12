"""Generalized isotropic elasticity with an independent shear modulus.

This constitutive form mirrors FEMaster's ``GENISO`` input and is kept separate
from ordinary isotropic elasticity because the additional shear modulus changes
the constitutive assumptions rather than being optional metadata.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .material_elasticity import Elasticity


class GeneralizedIsotropicElasticity(Elasticity):
    """Generalized isotropic elasticity defined by ``E``, ``nu`` and ``G``."""

    def __init__(
        self,
        youngs_modulus: float,
        poisson_ratio: float,
        shear_modulus: float,
    ) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio = float(poisson_ratio)
        self.shear_modulus = float(shear_modulus)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="GENISO"),
            csv((
                self.youngs_modulus,
                self.poisson_ratio,
                self.shear_modulus,
            )),
        ])
