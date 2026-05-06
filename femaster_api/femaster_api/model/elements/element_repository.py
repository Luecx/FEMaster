"""Repository for element objects."""

from __future__ import annotations

from dataclasses import replace
from typing import Iterator

from femaster_api.model.nodes import Node

from .element import Element
from .element_topology import TOPOLOGY_NODE_COUNTS


class ElementRepository:
    """Repository for finite elements in model order."""

    def __init__(self) -> None:
        self._elements: list[Element] = []

    def add(self, element: Element) -> Element:
        if not isinstance(element, Element):
            raise TypeError("element must be an Element")

        connectivity = tuple(_require_node(node) for node in element.nodes)
        valid_counts = TOPOLOGY_NODE_COUNTS[element.topology]
        if len(connectivity) not in valid_counts:
            expected = " or ".join(str(count) for count in valid_counts)
            raise ValueError(f"{element.topology.value} expects {expected} nodes, got {len(connectivity)}")

        item = replace(element, nodes=connectivity, id=len(self._elements))
        self._elements.append(item)
        return item

    def get(self, index: int) -> Element:
        try:
            return self._elements[index]
        except IndexError as exc:
            raise IndexError(f"unknown element index: {index}") from exc

    def __getitem__(self, index: int | slice) -> Element | tuple[Element, ...]:
        if isinstance(index, slice):
            return tuple(self._elements[index])
        return self.get(index)

    def has(self, element: Element) -> bool:
        return any(item is element for item in self._elements)

    def all(self) -> tuple[Element, ...]:
        return tuple(self._elements)

    def ids(self) -> tuple[int, ...]:
        return tuple(element.id for element in self._elements)

    def by_type(self):
        grouped = {}
        for element in self.all():
            grouped.setdefault(element.topology, []).append(element)
        return {key: tuple(value) for key, value in grouped.items()}

    def __iter__(self) -> Iterator[Element]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._elements)


def _require_node(value: Node) -> Node:
    if not isinstance(value, Node):
        raise TypeError(f"element connectivity must contain Node objects, got {type(value).__name__}")
    return value
