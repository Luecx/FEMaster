"""Direct generalized shell section with explicit object relationships.

The section stores a full 6x6 ABD matrix and 2x2 transverse-shear matrix in
row-major flattened order.  Its assignment target is an ``ElementRegion`` object;
optional material and orientation relationships are ``Material`` and
``CoordinateSystem`` objects.  Native names exist only at the serialization
boundary.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..material.material import Material
from ..region.region_element import ElementRegion
from .section import Section


class ABDShellSection(Section):
    """Direct generalized shell stiffness assignment."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        abd: Iterable[float],
        shear: Iterable[float],
        *,
        thickness: float = 1.0,
        material: Material | None = None,
        orientation: CoordinateSystem | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region)
        if material is not None and not isinstance(material, Material):
            raise TypeError("material must be a Material object")
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")

        self.abd = tuple(float(value) for value in abd)
        self.shear = tuple(float(value) for value in shear)
        self.thickness = float(thickness)
        self.material = material
        self.orientation = orientation
        self.csys_axis = int(csys_axis)

        if len(self.abd) != 36:
            raise ValueError("ABDShellSection requires exactly 36 ABD values")
        if len(self.shear) != 4:
            raise ValueError("ABDShellSection requires exactly 4 shear values")
        if self.csys_axis not in (1, 2, 3):
            raise ValueError("ABDShellSection csys_axis must be 1, 2 or 3")

    def export(self) -> str:
        values = (*self.abd, *self.shear)
        rows = [
            csv(values[start:start + 8])
            for start in range(0, len(values), 8)
        ]
        return block([
            keyword(
                "SHELLSECTION",
                TYPE="ABD",
                ELSET=self.element_region.name,
                MATERIAL=self.material.name if self.material is not None else None,
                THICKNESS=self.thickness,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
                CSYSAXIS=self.csys_axis,
            ),
            *rows,
        ])
