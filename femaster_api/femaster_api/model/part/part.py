"""Reusable Part-local model definition."""

from ..common.format import block, blocks, keyword
from ..common.named_object import NamedObject
from ..common.named_repository import NamedRepository
from ..mesh.element_repository import ElementRepository
from ..mesh.node_repository import NodeRepository
from ..mesh.surface import Surface
from ..region.region_repository import RegionRepository
from ..section.section_repository import SectionRepository


class Part(NamedObject):
    """Reusable Part owning local topology, regions, surfaces and sections."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.nodes = NodeRepository()
        self.elements = ElementRepository()
        self.regions = RegionRepository()
        self.surfaces = NamedRepository[Surface]()
        self.sections = SectionRepository()

    def export(self, *, root: bool = False) -> str:
        body = blocks((
            self.nodes.export(),
            self.elements.export(),
            self.regions.export(),
            "\n\n".join(surface.export() for surface in self.surfaces),
            self.sections.export(),
        ))
        if root:
            return body

        lines = [keyword("PART", NAME=self.name)]
        if body:
            lines.append(body)
        lines.append(keyword("ENDPART"))
        return block(lines)
