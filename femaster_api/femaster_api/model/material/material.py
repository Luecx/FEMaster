"""Named global material."""

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .elasticity import Elasticity


class Material(NamedObject):
    """Named global material composed from physical property definitions."""

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
