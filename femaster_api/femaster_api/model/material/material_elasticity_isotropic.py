"""Linear isotropic elasticity for a FEMaster material.

The model is defined by Young's modulus and Poisson ratio and exports the native
``*ELASTIC, TYPE=ISOTROPIC`` block.  Numerical values are stored exactly as
floating-point material parameters; unit consistency remains the caller's
responsibility as in the FEMaster input format.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .material_elasticity import Elasticity


class IsotropicElasticity(Elasticity):
    """Linear isotropic elasticity defined by ``E`` and ``nu``."""

    def __init__(self, youngs_modulus: float, poisson_ratio: float) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio = float(poisson_ratio)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="ISOTROPIC"),
            csv((self.youngs_modulus, self.poisson_ratio)),
        ])
