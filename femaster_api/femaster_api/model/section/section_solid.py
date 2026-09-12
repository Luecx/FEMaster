"""Continuum solid-section assignment.

A solid section links one part-local element region to one globally named
material and may additionally reference a coordinate-system orientation.  The
object maps directly to FEMaster's ``*SOLIDSECTION`` command and carries no
geometry-specific data row.
"""

from __future__ import annotations

from ..common.format import keyword
from .section_material import MaterialSection


class SolidSection(MaterialSection):
    """Material assignment for a solid-element region."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        orientation: str | None = None,
    ) -> None:
        super().__init__(name, element_region, material)
        self.orientation = orientation

    def export(self) -> str:
        return keyword(
            "SOLIDSECTION",
            ELSET=self.element_region,
            MATERIAL=self.material,
            ORIENTATION=self.orientation,
        )
