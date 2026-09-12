"""Container for the independent FEMaster region namespaces.

Node, element, surface and line regions intentionally use separate named
repositories because identical names may be meaningful in different domains.
The repository also exposes export phases because surface regions depend on
already-defined ``Surface`` objects: ``*NSET`` / ``*ELSET`` must be available
before topology that may reference them, while ``*SFSET`` must be written only
after the referenced ``*SURFACE`` definitions.

Line regions currently have no standalone native keyword and are therefore kept
in memory without inventing unsupported serialization syntax.
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

    def export_entity_regions(self) -> str:
        """Export node/element regions before dependent surface definitions."""

        rendered = [
            *(region.export() for region in self.nodes),
            *(region.export() for region in self.elements),
        ]
        return "\n\n".join(item for item in rendered if item)

    def export_surface_regions(self) -> str:
        """Export surface groups after their ``Surface`` objects are defined."""

        return "\n\n".join(
            item
            for item in (region.export() for region in self.surfaces)
            if item
        )

    def export(self) -> str:
        """Export all serializable regions when dependency phasing is irrelevant."""

        rendered = (
            self.export_entity_regions(),
            self.export_surface_regions(),
        )
        return "\n\n".join(item for item in rendered if item)
