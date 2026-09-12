"""Base section that references one global material definition.

Material-backed sections share the same semantic relation: a local element
region receives properties derived from a project-level material.  The class is
kept separate from ``Section`` so point-element properties can remain valid
without inventing a material reference they do not use.
"""

from __future__ import annotations

from .section import Section


class MaterialSection(Section):
    """Base class for section assignments referencing a material by name."""

    def __init__(self, name: str, element_region: str, material: str) -> None:
        super().__init__(name, element_region)
        self.material = str(material)
