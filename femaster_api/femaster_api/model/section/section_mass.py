"""Concentrated translational mass assigned to an ``ElementRegion`` object.

The property and topology remain separate, but the assignment itself is a direct
object relationship: ``MassSection.element_region`` is the concrete region, not
an ``ELSET`` string.  Export alone converts that relation to FEMaster's native
keyword representation.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..region.region_element import ElementRegion
from .section import Section


class MassSection(Section):
    """Isotropic concentrated mass property."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        mass: float,
    ) -> None:
        super().__init__(name, element_region)
        self.mass = float(mass)

    def export(self) -> str:
        return block([
            keyword("MASS", ELSET=self.element_region.name, TYPE="ISOTROPIC"),
            csv((self.mass,)),
        ])
