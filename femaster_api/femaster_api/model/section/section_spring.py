"""Ground-spring property assigned to ``SPRING1`` point elements.

The property stores one structural degree of freedom and one constant stiffness.
DOF validation is performed immediately because FEMaster structural point
springs address the six standard translational/rotational components.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section import Section


class SpringSection(Section):
    """Constant one-DOF ground spring property."""

    def __init__(
        self,
        name: str,
        element_region: str,
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
            keyword("SPRING", ELSET=self.element_region),
            csv((self.dof,)),
            csv((self.stiffness,)),
        ])
