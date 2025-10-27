from __future__ import annotations
from .section_base import Section
from ..materials.material import Material
from ..sets.elementset import ElementSet

class ShellSection(Section):
    def __init__(self, thickness: float, material: Material, elset: ElementSet):
        self.thickness = float(thickness)
        self.material = material
        self.elset = elset

    def to_femaster(self) -> str:
        return (
            f"*SHELL SECTION, ELSET={self.elset.name}, MATERIAL={self.material.name}\n"
            f"{self.thickness}"
        )

    def to_asami(self) -> str:
        raise NotImplementedError
