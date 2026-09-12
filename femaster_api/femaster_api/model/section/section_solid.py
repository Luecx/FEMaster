"""Three-dimensional continuum section assignment.

``SolidSection`` connects one element region to a global material and optionally
one material orientation.  It exports directly to FEMaster's native
``*SOLIDSECTION`` keyword and contains no element-specific serializer logic.
"""

from __future__ import annotations

from ..common.format import keyword
from .section_material import MaterialSection


class SolidSection(MaterialSection):
    """Material assignment for continuum solid elements."""

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
