"""Direct generalized shell section using ABD and transverse shear stiffness.

The object stores a full 6x6 ABD matrix and a 2x2 transverse shear matrix in
row-major flattened order.  Sizes are validated immediately because an
incomplete generalized section has no meaningful native representation.
Optional material and orientation references retain FEMaster's current
``SHELLSECTION TYPE=ABD`` capabilities.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .section import Section


class ABDShellSection(Section):
    """Direct generalized shell stiffness assignment."""

    def __init__(
        self,
        name: str,
        element_region: str,
        abd: Iterable[float],
        shear: Iterable[float],
        *,
        thickness: float = 1.0,
        material: str | None = None,
        orientation: str | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region)
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
                ELSET=self.element_region,
                MATERIAL=self.material,
                THICKNESS=self.thickness,
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            *rows,
        ])
