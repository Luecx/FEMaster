"""Base class for section/property assignments to concrete ``ElementRegion``s.

A section relates an element region to constitutive or concentrated properties.
The region is stored as the actual ``ElementRegion`` object; its semantic name is
used only when native ``ELSET=...`` syntax is emitted.  Repository positions and
string references never participate in the in-memory relationship.

Continuum, shell, beam and truss sections belong to a ``Part``.  Point-element
properties may additionally live in the assembly-level section repository where
FEMaster permits them.
"""

from __future__ import annotations

from ..common.named_object import NamedObject
from ..region.region_element import ElementRegion


class Section(NamedObject):
    """Base class for one named element-property assignment."""

    def __init__(self, name: str, element_region: ElementRegion) -> None:
        super().__init__(name)
        if not isinstance(element_region, ElementRegion):
            raise TypeError("element_region must be an ElementRegion object")
        self.element_region = element_region

    def export(self) -> str:
        raise NotImplementedError
