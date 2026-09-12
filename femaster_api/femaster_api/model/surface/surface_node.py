"""Node-based surface definition using ``Node`` / ``NodeRegion`` objects.

``NodeSurface`` represents native ``*SURFACE, TYPE=NODE``.  Membership is stored
as concrete node-domain objects instead of IDs or region-name strings.  Export
performs the only conversion back to native identifiers, keeping the editable
Python model strongly connected and type explicit.

Node surfaces are topology definitions in their own right and remain separate
from ``NodeRegion``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..node.node import Node
from ..region.region_node import NodeRegion
from .surface import Surface


class NodeSurface(Surface):
    """Named surface assembled from nodes or node regions."""

    def __init__(
        self,
        name: str,
        members: Iterable[Node | NodeRegion] = (),
    ) -> None:
        super().__init__(name)
        self.members: list[Node | NodeRegion] = []
        self.add(*members)

    def add(self, *members: Node | NodeRegion) -> "NodeSurface":
        """Append concrete node-domain objects and return this surface."""

        for member in members:
            if not isinstance(member, (Node, NodeRegion)):
                raise TypeError("NodeSurface members must be Node or NodeRegion")
        self.members.extend(members)
        return self

    def export(self) -> str:
        """Export one ``*SURFACE, TYPE=NODE`` block."""

        def value(target: Node | NodeRegion) -> int | str:
            return target.id if isinstance(target, Node) else target.name

        return block([
            keyword("SURFACE", NAME=self.name, TYPE="NODE"),
            *(csv((value(member),)) for member in self.members),
        ])
