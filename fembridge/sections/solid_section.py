from __future__ import annotations
from .section_base import Section
from ..materials.material import Material
from ..sets.elementset import ElementSet

class SolidSection(Section):
    def __init__(self, material: Material, elset: ElementSet):
        self.material = material
        self.elset = elset

    def to_femaster(self) -> str:
        return f"*SOLID SECTION, ELSET={self.elset.name}, MATERIAL={self.material.name}"

    def to_asami(self) -> str:
        raise NotImplementedError
