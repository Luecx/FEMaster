"""Isotropic point-mass property."""

from ..common.format import block, csv, keyword
from .section import Section


class MassSection(Section):
    """Isotropic translational mass assigned to MASS point elements."""

    def __init__(self, name: str, element_region: str, mass: float) -> None:
        super().__init__(name, element_region)
        self.mass = float(mass)

    def export(self) -> str:
        return block([
            keyword("MASS", ELSET=self.element_region, TYPE="ISOTROPIC"),
            csv((self.mass,)),
        ])
