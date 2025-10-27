# sections/beam_section.py
from __future__ import annotations
from typing import Tuple
from .section_base import Section
from ..materials.material import Material
from ..profiles.profile_base import Profile
from ..sets.elementset import ElementSet

class BeamSection(Section):
    def __init__(
            self,
            profile: Profile,
            material: Material,
            elset: ElementSet,
            n1: Tuple[float, float, float] = (1.0, 0.0, 0.0),
    ):
        self.profile = profile
        self.material = material
        self.elset = elset
        self.n1 = (float(n1[0]), float(n1[1]), float(n1[2]))

    def to_femaster(self) -> str:
        return "\n".join([
            f"*BEAM SECTION, ELSET={self.elset.name}, MATERIAL={self.material.name}, PROFILE={self.profile.name}",
            f"{self.n1[0]}, {self.n1[1]}, {self.n1[2]}",
        ])

    def to_asami(self) -> str:
        raise NotImplementedError
