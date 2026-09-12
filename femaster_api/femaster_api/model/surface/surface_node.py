"""Node-based surface definition.

``NodeSurface`` represents the native ``*SURFACE, TYPE=NODE`` form.  It is useful
when a FEMaster command expects surface semantics while the topology is supplied
as nodes or node-region references.  Membership remains ordered so export is
deterministic and imported decks can be reproduced without silent reordering.

Node surfaces are definitions in their own right and are therefore kept separate
from ``NodeRegion``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.typing import NodeReference
from .surface import Surface


class NodeSurface(Surface):
    """Named surface assembled from node or node-region references."""

    def __init__(self, name: str, members: Iterable[NodeReference] = ()) -> None:
        super().__init__(name)
        self.members: list[NodeReference] = list(members)

    def add(self, *members: NodeReference) -> "NodeSurface":
        """Append node references and return this surface."""

        self.members.extend(members)
        return self

    def export(self) -> str:
        """Export one ``*SURFACE, TYPE=NODE`` block."""

        return block([
            keyword("SURFACE", NAME=self.name, TYPE="NODE"),
            *(csv((member,)) for member in self.members),
        ])
