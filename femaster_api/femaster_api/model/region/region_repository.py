"""Container for independent FEMaster region namespaces."""

from ..common.named_repository import NamedRepository
from .element_region import ElementRegion
from .line_region import LineRegion
from .node_region import NodeRegion
from .region import Region
from .surface_region import SurfaceRegion


class RegionRepository:
    """Own node, element, surface and line region namespaces."""

    def __init__(self) -> None:
        self.nodes = NamedRepository[NodeRegion]()
        self.elements = NamedRepository[ElementRegion]()
        self.surfaces = NamedRepository[SurfaceRegion]()
        self.lines = NamedRepository[LineRegion]()

    def add(self, region: Region) -> Region:
        if isinstance(region, NodeRegion):
            return self.nodes.add(region)
        if isinstance(region, ElementRegion):
            return self.elements.add(region)
        if isinstance(region, SurfaceRegion):
            return self.surfaces.add(region)
        if isinstance(region, LineRegion):
            return self.lines.add(region)
        raise TypeError(f"unsupported region type: {type(region).__name__}")

    def export(self) -> str:
        rendered = [
            *(region.export() for region in self.nodes),
            *(region.export() for region in self.elements),
        ]
        return "\n\n".join(item for item in rendered if item)
