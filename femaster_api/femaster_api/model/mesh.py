"""Part-local mesh entities and repositories.

Nodes and elements preserve the sparse integer identifiers supplied by the user.
Repositories own deterministic model order but never renumber topology. Concrete
element classes carry the FEMaster TYPE token and validate connectivity size.
"""

from __future__ import annotations

from .typing import ElementReference
from .._format import block, csv, keyword
from ..repository import IdRepository


class Node:
    """Part-local finite-element node."""

    def __init__(self, id: int, x: float, y: float, z: float = 0.0) -> None:
        self.id = int(id)
        self.x  = float(x)
        self.y  = float(y)
        self.z  = float(z)

    def to_femaster(self) -> str:
        """Return this node as one FEMaster NODE data row."""

        return csv((self.id, self.x, self.y, self.z))


class NodeRepository(IdRepository[Node]):
    """Repository of part-local nodes keyed by FEMaster node id."""

    def to_femaster(self) -> str:
        """Return all nodes as one FEMaster NODE block."""

        if not self._items:
            return ""
        return block([keyword("NODE"), *(node.to_femaster() for node in self)])


class Element:
    """Base class for part-local finite elements with fixed connectivity."""

    type_name = "ELEMENT"
    node_count: int | None = None

    def __init__(self, id: int, nodes: tuple[int, ...] | list[int]) -> None:
        self.id    = int(id)
        self.nodes = tuple(int(node) for node in nodes)

        if self.node_count is not None and len(self.nodes) != self.node_count:
            raise ValueError(
                f"{self.type_name} requires {self.node_count} nodes, got {len(self.nodes)}"
            )

    def to_femaster(self) -> str:
        """Return this element as one FEMaster ELEMENT connectivity row."""

        return csv((self.id, *self.nodes))


class C3D4(Element):
    type_name, node_count = "C3D4", 4


class C3D5(Element):
    type_name, node_count = "C3D5", 5


class C3D6(Element):
    type_name, node_count = "C3D6", 6


class C3D8(Element):
    type_name, node_count = "C3D8", 8


class C3D8R(Element):
    type_name, node_count = "C3D8R", 8


class C3D10(Element):
    type_name, node_count = "C3D10", 10


class C3D15(Element):
    type_name, node_count = "C3D15", 15


class C3D20(Element):
    type_name, node_count = "C3D20", 20


class C3D20R(Element):
    type_name, node_count = "C3D20R", 20


class C2D3(Element):
    type_name, node_count = "C2D3", 3


class C2D4(Element):
    type_name, node_count = "C2D4", 4


class C2D6(Element):
    type_name, node_count = "C2D6", 6


class C2D8(Element):
    type_name, node_count = "C2D8", 8


class S3(Element):
    type_name, node_count = "S3", 3


class S4(Element):
    type_name, node_count = "S4", 4


class S6(Element):
    type_name, node_count = "S6", 6


class S8(Element):
    type_name, node_count = "S8", 8


class B33(Element):
    type_name, node_count = "B33", 2


class T3D2(Element):
    type_name, node_count = "T3D2", 2


ELEMENT_TYPES: dict[str, type[Element]] = {
    item.type_name: item
    for item in (
        C3D4, C3D5, C3D6, C3D8, C3D8R, C3D10, C3D15, C3D20, C3D20R,
        C2D3, C2D4, C2D6, C2D8,
        S3, S4, S6, S8,
        B33, T3D2,
    )
}


class ElementRepository(IdRepository[Element]):
    """Repository of part-local elements keyed by FEMaster element id."""

    def to_femaster(self) -> str:
        """Return elements grouped into deterministic FEMaster TYPE blocks."""

        groups: dict[str, list[Element]] = {}
        for element in self:
            groups.setdefault(element.type_name, []).append(element)

        result: list[str] = []
        for type_name, elements in groups.items():
            result.append(
                block([
                    keyword("ELEMENT", TYPE=type_name),
                    *(element.to_femaster() for element in elements),
                ])
            )
        return "\n\n".join(result)


class Surface:
    """Named surface assembled from element-side references."""

    def __init__(self, name: str) -> None:
        from ..repository import NamedObject

        self._named = NamedObject(name)
        self.entries: list[tuple[ElementReference, int]] = []

    @property
    def name(self) -> str:
        return self._named.name

    def add(self, element: ElementReference, side: int) -> "Surface":
        """Append one element or ELSET side reference."""

        self.entries.append((element, int(side)))
        return self

    def to_femaster(self) -> str:
        """Return this surface as one FEMaster SURFACE block."""

        return block([
            keyword("SURFACE", NAME=self.name),
            *(csv((element, side)) for element, side in self.entries),
        ])
