"""Axial truss section relating concrete region and material objects.

``TrussSection`` stores an ``ElementRegion`` and ``Material`` directly together
with a strictly positive cross-sectional area.  Only native export converts the
two object relationships to their semantic names.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..material.material import Material
from ..region.region_element import ElementRegion
from .section_material import MaterialSection


class TrussSection(MaterialSection):
    """Material-backed truss section defined by cross-sectional area."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        material: Material,
        area: float,
    ) -> None:
        super().__init__(name, element_region, material)
        self.area = float(area)
        if self.area <= 0.0:
            raise ValueError("TrussSection area must be positive")

    def export(self) -> str:
        return block([
            keyword(
                "TRUSSSECTION",
                ELSET=self.element_region.name,
                MATERIAL=self.material.name,
            ),
            csv((self.area,)),
        ])
