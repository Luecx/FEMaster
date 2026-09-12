"""Axial truss section assignment.

``TrussSection`` associates an element region with one global material and a
strictly positive cross-sectional area.  The object maps directly to
``*TRUSSSECTION`` and performs the only local physical validation required for
that native definition.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section_material import MaterialSection


class TrussSection(MaterialSection):
    """Material-backed truss section defined by cross-sectional area."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
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
                ELSET=self.element_region,
                MATERIAL=self.material,
            ),
            csv((self.area,)),
        ])
