"""Ordered repository for model nodes."""

from __future__ import annotations

from typing import Iterator

from .node import Node


class NodeRepository:
    """Container of nodes in model order."""

    def __init__(self) -> None:
        self._nodes: list[Node] = []

    def add(self, node: Node) -> Node:
        """Add a node and return the stored object."""

        if not isinstance(node, Node):
            raise TypeError("node must be a Node")
        item = Node(float(node.x), float(node.y), float(node.z))
        self._nodes.append(item)
        return item

    def get(self, index: int) -> Node:
        try:
            return self._nodes[index]
        except IndexError as exc:
            raise IndexError(f"unknown node index: {index}") from exc

    def __getitem__(self, index: int | slice) -> Node | tuple[Node, ...]:
        if isinstance(index, slice):
            return tuple(self._nodes[index])
        return self.get(index)

    def has(self, node: Node) -> bool:
        return any(item is node for item in self._nodes)

    def __contains__(self, node: object) -> bool:
        return isinstance(node, Node) and self.has(node)

    def all(self) -> tuple[Node, ...]:
        return tuple(self._nodes)

    def __iter__(self) -> Iterator[Node]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._nodes)
