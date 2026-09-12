"""Base representation of one part-local finite element.

Every concrete FEMaster element owns a sparse element ID and a connectivity made
of actual ``Node`` objects.  Node IDs are therefore never used as an in-memory
cross-object reference: they are read only when parsing a deck and written only
when serializing the element row.  This makes connectivity navigable and prevents
a Python model from containing dangling integer references after construction.

Concrete element classes define only their native ``TYPE`` token and required
node count.  Grouping several element rows under one ``*ELEMENT`` header remains
the responsibility of ``ElementRepository``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import csv
from ..node.node import Node


class Element:
    """Base finite element with persistent ID and object-valued connectivity."""

    type_name = "ELEMENT"
    node_count: int | None = None

    def __init__(self, id: int, nodes: Iterable[Node]) -> None:
        self.id = int(id)
        self.nodes = tuple(nodes)

        if not all(isinstance(node, Node) for node in self.nodes):
            raise TypeError("element connectivity must contain Node objects")

        if self.node_count is not None and len(self.nodes) != self.node_count:
            raise ValueError(
                f"{self.type_name} requires {self.node_count} nodes, "
                f"got {len(self.nodes)}"
            )

    def export(self) -> str:
        """Export this element using the persistent IDs of its connected nodes."""

        return csv((self.id, *(node.id for node in self.nodes)))
