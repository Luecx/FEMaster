"""Generalized isotropic elasticity."""

from ..common.format import block, csv, keyword
from .elasticity import Elasticity


class GeneralizedIsotropicElasticity(Elasticity):
    """Isotropic elasticity with independent shear modulus."""

    def __init__(self, youngs_modulus: float, poisson_ratio: float, shear_modulus: float) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio = float(poisson_ratio)
        self.shear_modulus = float(shear_modulus)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="GENISO"),
            csv((self.youngs_modulus, self.poisson_ratio, self.shear_modulus)),
        ])
