"""Container for the independent FEMaster region namespaces.

Node, element, surface and line regions intentionally use separate named
repositories because identical names may be meaningful in different entity
domains and FEMaster itself resolves them through domain-specific contexts.
``RegionRepository`` is therefore a small facade rather than one flattened set
table.

Only domains with a direct native keyword are emitted by ``export``.  A line
region can still be represented in memory without inventing unsupported syntax.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .region import Region
from .region_element import ElementRegion
from .region_line import LineRegion
from .region_node import NodeRegion
from .region_surface import SurfaceRegion


class RegionRepository:
    """Own independent named repositories for every supported region domain."""

    def __init__(self) -> None:
        self.nodes = NamedRepository[NodeRegion]()
        self.elements = NamedRepository[ElementRegion]()
        self.surfaces = NamedRepository[SurfaceRegion]()
        self.lines = NamedRepository[LineRegion]()

    def add(self, region: Region) -> Region:
        """Insert a region into the repository matching its concrete domain."""

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
        """Export all region domains with a native standalone representation."""

        rendered = [
            *(region.export() for region in self.nodes),
            *(region.export() for region in self.elements),
            *(region.export() for region in self.surfaces),
        ]
        return "\n\n".join(item for item in rendered if item)
