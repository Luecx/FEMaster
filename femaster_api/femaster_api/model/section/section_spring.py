"""Ground-spring property assigned to a concrete ``ElementRegion`` object.

The property stores one structural degree of freedom and one constant stiffness.
Its assignment target is a real region object; the region name is produced only
for native ``*SPRING`` serialization.  DOF validation is performed immediately.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..region.region_element import ElementRegion
from .section import Section


class SpringSection(Section):
    """Constant one-DOF ground spring property."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        dof: int,
        stiffness: float,
    ) -> None:
        super().__init__(name, element_region)
        self.dof = int(dof)
        self.stiffness = float(stiffness)
        if self.dof < 1 or self.dof > 6:
            raise ValueError("SpringSection dof must be between 1 and 6")

    def export(self) -> str:
        return block([
            keyword("SPRING", ELSET=self.element_region.name),
            csv((self.dof,)),
            csv((self.stiffness,)),
        ])
