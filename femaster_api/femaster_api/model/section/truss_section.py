"""Truss section assignment."""

from ..common.format import block, csv, keyword
from .material_section import MaterialSection


class TrussSection(MaterialSection):
    """Axial truss section defined by cross-sectional area."""

    def __init__(self, name: str, element_region: str, material: str, area: float) -> None:
        super().__init__(name, element_region, material)
        self.area = float(area)
        if self.area <= 0.0:
            raise ValueError("TrussSection area must be positive")

    def export(self) -> str:
        return block([
            keyword("TRUSSSECTION", ELSET=self.element_region, MATERIAL=self.material),
            csv((self.area,)),
        ])
