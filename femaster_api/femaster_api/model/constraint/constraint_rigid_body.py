"""Rigid-body-mode constraint/removal on one concrete ``ElementRegion``.

The constraint stores the actual element-region object affected by ``*RBM``.
The native ``ELSET`` name is derived from that object only during serialization,
which prevents unresolved region-name strings from entering the Python model.
"""

from __future__ import annotations

from ..common.format import keyword
from ..region.region_element import ElementRegion
from .constraint import Constraint


class RigidBodyConstraint(Constraint):
    """Rigid-body constraint acting on one element region."""

    def __init__(self, element_region: ElementRegion) -> None:
        if not isinstance(element_region, ElementRegion):
            raise TypeError("element_region must be an ElementRegion object")
        self.element_region = element_region

    def export(self) -> str:
        return keyword("RBM", ELSET=self.element_region.name)
