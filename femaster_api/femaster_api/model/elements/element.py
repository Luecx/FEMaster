"""Element data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.nodes import Node

from .element_topology import ElementTopology, FEMASTER_ELEMENT_TYPES


@dataclass(frozen=True, slots=True, eq=False)
class Element:
    """Finite element connectivity storing node objects."""

    topology: ElementTopology
    nodes: tuple[Node, ...]
    id: int | None = None

    @property
    def femaster_type(self) -> str:
        return FEMASTER_ELEMENT_TYPES[self.topology]

    @property
    def node_ids(self) -> tuple[int, ...]:
        return tuple(_node_id(node) for node in self.nodes)


def _node_id(node: Node) -> int:
    if node.id is None:
        raise ValueError("node has no repository id")
    return node.id
