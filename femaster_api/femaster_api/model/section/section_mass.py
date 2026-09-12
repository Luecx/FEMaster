"""Concentrated translational mass assigned to ``MASS`` point elements.

The topology and property remain separate: ``MassElement`` supplies the one-node
element while ``MassSection`` assigns the scalar mass to an element region.  The
same property keyword may be legal in part or assembly scope; repository
ownership determines where it is exported.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section import Section


class MassSection(Section):
    """Isotropic concentrated mass property."""

    def __init__(self, name: str, element_region: str, mass: float) -> None:
        super().__init__(name, element_region)
        self.mass = float(mass)

    def export(self) -> str:
        return block([
            keyword("MASS", ELSET=self.element_region, TYPE="ISOTROPIC"),
            csv((self.mass,)),
        ])
