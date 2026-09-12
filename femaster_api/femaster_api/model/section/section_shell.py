"""Integrated shell section with constant thickness.

This class represents the material-integrated shell form, not a direct ABD
stiffness.  Thickness, material orientation and coordinate-system axis are
stored explicitly so the object corresponds one-to-one with the native
``*SHELLSECTION, TYPE=INTEGRATED`` definition.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from .section_material import MaterialSection


class ShellSection(MaterialSection):
    """Material-integrated shell section with constant thickness."""

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
        if self.csys_axis not in (1, 2, 3):
            raise ValueError("ShellSection csys_axis must be 1, 2 or 3")

    def export(self) -> str:
        return block([
            keyword(
                "SHELLSECTION",
                TYPE="INTEGRATED",
                ELSET=self.element_region,
                MATERIAL=self.material,
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            csv((self.thickness,)),
        ])
