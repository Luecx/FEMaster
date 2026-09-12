"""Generalized shell section defined directly by ABD and transverse-shear terms.

FEMaster's ABD shell-section form consumes exactly 36 coupling/stiffness values
followed by four transverse-shear coefficients.  The Python class stores those
groups separately for readability, validates their dimensions immediately and
emits the native forty-value sequence unchanged.

Material and orientation names remain optional because generalized ABD input can
be fully constitutive on its own while still supporting semantic metadata needed
by downstream shell calculations.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .section import Section


class ABDShellSection(Section):
    """Generalized 6x6 ABD plus 2x2 transverse-shear shell section."""

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
        rows = [csv(values[start:start + 8]) for start in range(0, len(values), 8)]
        return block([
            keyword(
                "SHELLSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                TYPE="ABD",
                THICKNESS=self.thickness,
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            *rows,
        ])
