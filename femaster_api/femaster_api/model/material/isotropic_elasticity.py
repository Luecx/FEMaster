"""Linear isotropic elasticity."""

from ..common.format import block, csv, keyword
from .elasticity import Elasticity


class IsotropicElasticity(Elasticity):
    """Young's modulus and Poisson-ratio elasticity."""

    def __init__(self, youngs_modulus: float, poisson_ratio: float) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio = float(poisson_ratio)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="ISOTROPIC"),
            csv((self.youngs_modulus, self.poisson_ratio)),
        ])
