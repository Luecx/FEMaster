"""Named global material definition.

``Material`` is a project-level definition shared by any number of part-local
sections.  Physical properties remain explicit attributes and constituent
objects rather than being hidden in a dictionary or serializer registry.
Currently the class composes elasticity, density and scalar thermal expansion in
the same dependency order expected by FEMaster.

The material name is immutable because sections reference it semantically.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .material_elasticity import Elasticity


class Material(NamedObject):
    """Named material composed from independent constitutive properties."""

    def __init__(
        self,
        name: str,
        *,
        elasticity: Elasticity | None = None,
        density: float | None = None,
        thermal_expansion: float | None = None,
    ) -> None:
        super().__init__(name)
        self.elasticity = elasticity
        self.density = None if density is None else float(density)
        self.thermal_expansion = (
            None if thermal_expansion is None else float(thermal_expansion)
        )

    def export(self) -> str:
        """Export the complete native ``*MATERIAL`` definition."""

        lines = [keyword("MATERIAL", NAME=self.name)]

        if self.elasticity is not None:
            lines.append(self.elasticity.export())

        if self.density is not None:
            lines.extend([keyword("DENSITY"), csv((self.density,))])

        if self.thermal_expansion is not None:
            lines.extend([
                keyword("THERMALEXPANSION"),
                csv((self.thermal_expansion,)),
            ])

        return block(lines)
