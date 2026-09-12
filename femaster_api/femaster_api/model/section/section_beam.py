"""Beam section assignment referencing material, profile and local direction.

The section keeps the beam profile and material as global semantic names while
the orientation vector defines the local section direction required by the
native ``*BEAMSECTION`` data row.  A zero or malformed direction is rejected at
construction time rather than producing an invalid deck.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .section_material import MaterialSection


class BeamSection(MaterialSection):
    """Beam section with named profile and explicit local orientation."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        profile: str,
        orientation: Iterable[float],
    ) -> None:
        super().__init__(name, element_region, material)
        self.profile = str(profile)
        self.orientation = tuple(float(value) for value in orientation)

        if len(self.orientation) != 3:
            raise ValueError("BeamSection orientation requires exactly 3 values")
        if all(value == 0.0 for value in self.orientation):
            raise ValueError("BeamSection orientation must be non-zero")

    def export(self) -> str:
        return block([
            keyword(
                "BEAMSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                PROFILE=self.profile,
            ),
            csv(self.orientation),
        ])
