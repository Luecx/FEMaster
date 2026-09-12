"""Ground-spring point property."""

from ..common.format import block, csv, keyword
from .section import Section


class SpringSection(Section):
    """Constant ground stiffness assigned to SPRING1 point elements."""

    def __init__(self, name: str, element_region: str, dof: int, stiffness: float) -> None:
        super().__init__(name, element_region)
        self.dof = int(dof)
        self.stiffness = float(stiffness)
        if not 1 <= self.dof <= 6:
            raise ValueError("SpringSection dof must be between 1 and 6")

    def export(self) -> str:
        return block([
            keyword("SPRING", ELSET=self.element_region),
            csv((self.dof,)),
            csv((self.stiffness,)),
        ])
