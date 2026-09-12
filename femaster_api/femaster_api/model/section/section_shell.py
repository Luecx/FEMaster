"""Integrated shell section with concrete region/material/orientation objects.

This class represents the material-integrated shell form, not a direct ABD
stiffness.  ``ElementRegion`` and ``Material`` are mandatory object references;
an optional material orientation is a ``CoordinateSystem`` object.  Thickness
and coordinate-system axis remain intrinsic numerical section data.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..material.material import Material
from ..region.region_element import ElementRegion
from .section_material import MaterialSection


class ShellSection(MaterialSection):
    """Material-integrated shell section with constant thickness."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        material: Material,
        thickness: float,
        orientation: CoordinateSystem | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region, material)
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")
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
                ELSET=self.element_region.name,
                MATERIAL=self.material.name,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
                CSYSAXIS=self.csys_axis,
            ),
            csv((self.thickness,)),
        ])
