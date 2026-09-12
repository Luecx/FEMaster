"""Beam section assignment."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .material_section import MaterialSection


class BeamSection(MaterialSection):
    """Beam section referencing a global profile and n1 direction."""

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
