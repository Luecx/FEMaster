"""Base class for section and point-property assignments.

A section associates an element region with constitutive or concentrated
properties.  The Python object has its own immutable name for repository lookup
although several FEMaster section keywords identify the assignment only by
``ELSET``.  Cross-object references therefore remain explicit strings rather
than repository positions.

Continuum, shell, beam and truss sections belong to a ``Part``.  Only the
point-element properties explicitly allowed by FEMaster may additionally live in
the assembly-level section repository on ``Project``.
"""

from __future__ import annotations

from ..common.named_object import NamedObject


class Section(NamedObject):
    """Base class for one named element-property assignment."""

    def __init__(self, name: str, element_region: str) -> None:
        super().__init__(name)
        self.element_region = str(element_region)

    def export(self) -> str:
        raise NotImplementedError
