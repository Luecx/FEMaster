"""Base section that relates an ``ElementRegion`` to a concrete ``Material``.

Material-backed sections share one semantic relationship: an element region gets
properties from one project-level material object.  Both relationships are held
as real objects and only reduced to ``ELSET`` / ``MATERIAL`` names during native
serialization.

Point-element properties remain derived directly from ``Section`` because they
do not require a material.
"""

from __future__ import annotations

from ..material.material import Material
from ..region.region_element import ElementRegion
from .section import Section


class MaterialSection(Section):
    """Base class for section assignments backed by one material object."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        material: Material,
    ) -> None:
        super().__init__(name, element_region)
        if not isinstance(material, Material):
            raise TypeError("material must be a Material object")
        self.material = material
