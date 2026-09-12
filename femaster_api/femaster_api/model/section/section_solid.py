"""Three-dimensional continuum section with object-valued model relations.

``SolidSection`` connects an ``ElementRegion`` to a ``Material`` and optionally a
``CoordinateSystem``.  Those objects are retained directly; FEMaster's semantic
names are derived only when ``*SOLIDSECTION`` is exported.
"""

from __future__ import annotations

from ..common.format import keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..material.material import Material
from ..region.region_element import ElementRegion
from .section_material import MaterialSection


class SolidSection(MaterialSection):
    """Material assignment for continuum solid elements."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        material: Material,
        orientation: CoordinateSystem | None = None,
    ) -> None:
        super().__init__(name, element_region, material)
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")
        self.orientation = orientation

    def export(self) -> str:
        return keyword(
            "SOLIDSECTION",
            ELSET=self.element_region.name,
            MATERIAL=self.material.name,
            ORIENTATION=(
                self.orientation.name if self.orientation is not None else None
            ),
        )
