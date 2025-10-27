
from __future__ import annotations
from typing import Optional, List

from .elasticity_iso import ElasticityIsotropic
from .elasticity_ortho import ElasticityOrthotropic
from .density import Density
from .expansion_iso import ThermalExpansionIsotropic
from .conductivity_iso import ConductivityIsotropic
from .specific_heat import SpecificHeat


class Material:
    """
    A single material definition composed from property classes.
    """

    def __init__(
        self,
        name: str,
        *,
        elasticity: ElasticityIsotropic | ElasticityOrthotropic | None = None,
        density: Density | None = None,
        expansion: ThermalExpansionIsotropic | None = None,
        conductivity: ConductivityIsotropic | None = None,
        specific_heat: SpecificHeat | None = None,
    ):
        if not name:
            raise ValueError("Material name must be non-empty.")
        self.name = str(name)

        self.elasticity = elasticity
        self.density = density
        self.expansion = expansion
        self.conductivity = conductivity
        self.specific_heat = specific_heat

    def set_elasticity(self, e): self.elasticity = e; return self
    def set_density(self, d): self.density = d; return self
    def set_expansion(self, a): self.expansion = a; return self
    def set_conductivity(self, k): self.conductivity = k; return self
    def set_specific_heat(self, cp): self.specific_heat = cp; return self

    def to_femaster(self) -> str:
        lines: List[str] = [f"*MATERIAL, NAME={self.name}"]
        if self.elasticity: lines.extend(self.elasticity.to_femaster())
        if self.density: lines.extend(self.density.to_femaster())
        if self.expansion: lines.extend(self.expansion.to_femaster())
        if self.conductivity: lines.extend(self.conductivity.to_femaster())
        if self.specific_heat: lines.extend(self.specific_heat.to_femaster())
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError("Material.to_asami not implemented yet.")
