"""Truss-section assignment with scalar cross-sectional area.

The section links a truss element region to a material and one positive area.
The value is validated before export so malformed zero/negative cross sections
cannot reach the native ``*TRUSSSECTION`` block.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section_material import MaterialSection


class TrussSection(MaterialSection):
    """Material and area assignment for a truss-element region."""

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
