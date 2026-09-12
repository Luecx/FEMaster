"""Typed named FEMaster regions.

Node and element regions have direct NSET/ELSET input representations. Surface
and line regions are retained as semantic region types, but their topology must
be constructed through SURFACE definitions; they therefore are not serialized
as invented SURFACESET/LINESET keywords.
"""

from __future__ import annotations

from collections.abc import Iterable

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository
from .typing import EntityReference


class Region(NamedObject):
    """Base class for one named FEMaster entity region."""

    keyword_name: str | None = None

    def __init__(self, name: str, members: Iterable[EntityReference] = ()) -> None:
        super().__init__(name)
        self.members: list[EntityReference] = list(members)

    def add(self, *members: EntityReference) -> "Region":
        """Append entity references and return this region."""

        self.members.extend(members)
        return self

    def export(self) -> str:
        """Export this region when a direct set keyword exists."""

        if self.keyword_name is None:
            raise NotImplementedError(
                f"{type(self).__name__} requires a topology definition and has no standalone keyword"
            )

        rows: list[str] = []
        for start in range(0, len(self.members), 16):
            rows.append(csv(self.members[start:start + 16]))

        return block([keyword(self.keyword_name, NAME=self.name), *rows])


class NodeRegion(Region):
    """Named node region written through NSET."""

    keyword_name = "NSET"


class ElementRegion(Region):
    """Named element region written through ELSET."""

    keyword_name = "ELSET"


class SurfaceRegion(Region):
    """Semantic region of materialized surface identifiers.

    Surface topology is exported through ``Surface`` definitions rather than a
    standalone set keyword.
    """


class LineRegion(Region):
    """Semantic region of materialized line identifiers.

    Line topology is produced by element-boundary SURFACE definitions and has no
    independent line-set input keyword.
    """


class RegionRepository:
    """Container for the four independent FEMaster region namespaces."""

    def __init__(self) -> None:
        self.nodes    = NamedRepository[NodeRegion]()
        self.elements = NamedRepository[ElementRegion]()
        self.surfaces = NamedRepository[SurfaceRegion]()
        self.lines    = NamedRepository[LineRegion]()

    def add(self, region: Region) -> Region:
        """Insert a typed region into its matching namespace."""

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
        """Export directly representable NSET and ELSET definitions."""

        rendered = [
            *(region.export() for region in self.nodes),
            *(region.export() for region in self.elements),
        ]
        return "\n\n".join(item for item in rendered if item)
