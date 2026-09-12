"""Reusable FEMaster part definition.

A ``Part`` is the sole owner of part-local nodes, elements, regions, surfaces and
section assignments.  This mirrors the C++ model before compilation and keeps
local IDs meaningful across repeated instances.  Shared definitions such as
materials, profiles and amplitudes remain on ``Project`` and are referenced by
semantic name.

The same class also represents the implicit root/default part.  The owning
``PartRepository`` decides whether export wraps the content in ``*PART`` /
``*ENDPART`` or writes it directly in root scope.
"""

from __future__ import annotations

from ..common.format import block, blocks, keyword
from ..common.named_object import NamedObject
from ..element.element_repository import ElementRepository
from ..node.node_repository import NodeRepository
from ..region.region_repository import RegionRepository
from ..section.section_repository import SectionRepository
from ..surface.surface_repository import SurfaceRepository


class Part(NamedObject):
    """Reusable part-local finite-element model definition."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.nodes = NodeRepository()
        self.elements = ElementRepository()
        self.regions = RegionRepository()
        self.surfaces = SurfaceRepository()
        self.sections = SectionRepository()

    def export(self, *, root: bool = False) -> str:
        """Export this part either in root scope or as an explicit ``*PART``."""

        body = blocks((
            self.nodes.export(),
            self.elements.export(),
            self.regions.export(),
            self.surfaces.export(),
            self.sections.export(),
        ))

        if root:
            return body

        lines = [keyword("PART", NAME=self.name)]
        if body:
            lines.append(body)
        lines.append(keyword("ENDPART"))
        return block(lines)
