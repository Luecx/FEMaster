"""Integrated shell section with material, thickness and orientation controls.

The section represents FEMaster's ordinary integrated shell formulation.  It
references one global material, stores one positive thickness and may reference
a named orientation.  ``csys_axis`` selects which local coordinate-system axis
is used and is validated against FEMaster's supported values 1..3.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section_material import MaterialSection


class ShellSection(MaterialSection):
    """Integrated shell section with one material and one thickness."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        thickness: float,
        orientation: str | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region, material)
        self.thickness = float(thickness)
        self.orientation = orientation
        self.csys_axis = int(csys_axis)

        if self.thickness <= 0.0:
            raise ValueError("ShellSection thickness must be positive")
        if self.csys_axis not in (1, 2, 3):
            raise ValueError("ShellSection csys_axis must be 1, 2 or 3")

    def export(self) -> str:
        return block([
            keyword(
                "SHELLSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                TYPE="INTEGRATED",
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            csv((self.thickness,)),
        ])
