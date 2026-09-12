"""Named global beam-profile definition.

A profile stores the scalar section properties consumed by FEMaster beam
sections.  It is independent of material and element-region assignment, allowing
one profile to be reused by multiple ``BeamSection`` objects.  Values are kept
in the native profile order so export remains transparent and reversible.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject


class Profile(NamedObject):
    """General beam profile expressed by FEMaster's scalar section properties."""

    def __init__(
        self,
        name: str,
        area: float,
        iy: float,
        iz: float,
        j: float,
        iyz: float = 0.0,
        ey: float = 0.0,
        ez: float = 0.0,
        refy: float = 0.0,
        refz: float = 0.0,
    ) -> None:
        super().__init__(name)
        self.values = tuple(
            float(value)
            for value in (area, iy, iz, j, iyz, ey, ez, refy, refz)
        )

    def export(self) -> str:
        return block([
            keyword("PROFILE", NAME=self.name),
            csv(self.values),
        ])
