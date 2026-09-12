"""Solid section assignment."""

from ..common.format import keyword
from .material_section import MaterialSection


class SolidSection(MaterialSection):
    """Three-dimensional continuum section assignment."""

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
