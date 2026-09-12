"""Beam section relating region, material and profile objects explicitly.

``BeamSection`` keeps the concrete ``ElementRegion``, ``Material`` and ``Profile``
objects used by the section.  The local orientation vector remains numerical
section data rather than a model-object reference.  Native names are generated
only while serializing ``*BEAMSECTION``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..material.material import Material
from ..profile.profile import Profile
from ..region.region_element import ElementRegion
from .section_material import MaterialSection


class BeamSection(MaterialSection):
    """Beam section with concrete profile and explicit local direction."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        material: Material,
        profile: Profile,
        orientation: Iterable[float],
    ) -> None:
        super().__init__(name, element_region, material)
        if not isinstance(profile, Profile):
            raise TypeError("profile must be a Profile object")
        self.profile = profile
        self.orientation = tuple(float(value) for value in orientation)

        if len(self.orientation) != 3:
            raise ValueError("BeamSection orientation requires exactly 3 values")
        if all(value == 0.0 for value in self.orientation):
            raise ValueError("BeamSection orientation must be non-zero")

    def export(self) -> str:
        return block([
            keyword(
                "BEAMSECTION",
                ELSET=self.element_region.name,
                MATERIAL=self.material.name,
                PROFILE=self.profile.name,
            ),
            csv(self.orientation),
        ])
