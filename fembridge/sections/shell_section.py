from __future__ import annotations
from .section_base import Section
from ..materials.material import Material
from ..sets.elementset import ElementSet

class ShellSection(Section):
    def __init__(self, thickness: float, material: Material, elset: ElementSet, orientation=None):
        self.thickness = float(thickness)
        self.material = material
        self.elset = elset
        self.orientation = orientation

    def to_femaster(self) -> str:
        orientation = ""
        if self.orientation is not None:
            orientation_name = getattr(self.orientation, "name", self.orientation)
            orientation = f", ORIENTATION={orientation_name}"
        return (
            f"*SHELL SECTION, ELSET={self.elset.name}, MATERIAL={self.material.name}{orientation}\n"
            f"{self.thickness}"
        )

    def to_asami(self) -> str:
        raise NotImplementedError
