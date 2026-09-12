"""Part-local mesh entities and repositories.

Nodes and elements preserve the sparse integer identifiers supplied by the user.
Repositories own deterministic model order but never renumber topology. Concrete
element classes mirror the TYPE names registered by the current FEMaster parser
and validate connectivity size before export.
"""

from __future__ import annotations

from .._format import block, csv, keyword
from ..repository import IdRepository, NamedObject
from .typing import ElementReference


class Node:
    """Part-local finite-element node."""

    def __init__(self, id: int, x: float, y: float, z: float = 0.0) -> None:
        self.id = int(id)
        self.x  = float(x)
        self.y  = float(y)
        self.z  = float(z)

    def export(self) -> str:
        """Export this node as one NODE data row."""

        return csv((self.id, self.x, self.y, self.z))


class NodeRepository(IdRepository[Node]):
    """Repository of part-local nodes keyed by FEMaster node id."""

    def export(self) -> str:
        """Export all nodes as one NODE block."""

        if not self._items:
            return ""
        return block([keyword("NODE"), *(node.export() for node in self)])


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

    def export(self) -> str:
        """Export this element as one ELEMENT connectivity row."""

        return csv((self.id, *self.nodes))


# -----------------------------------------------------------------------------
# Solid elements
# -----------------------------------------------------------------------------


class C3D4(Element):
    type_name, node_count = "C3D4", 4


class C3D5(Element):
    """Five-node pyramid accepted by FEMaster and expanded internally to C3D8."""

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


# -----------------------------------------------------------------------------
# Beam and truss elements
# -----------------------------------------------------------------------------


class B33(Element):
    type_name, node_count = "B33", 2


class T3(Element):
    """Two-node truss using the native short FEMaster TYPE name."""

    type_name, node_count = "T3", 2


class T3D2(Element):
    """Abaqus-compatible alias for the same two-node truss formulation."""

    type_name, node_count = "T3D2", 2


# -----------------------------------------------------------------------------
# Shell elements
# -----------------------------------------------------------------------------


class S3(Element):
    type_name, node_count = "S3", 3


class S4(Element):
    type_name, node_count = "S4", 4


class S6(Element):
    type_name, node_count = "S6", 6


class S8(Element):
    type_name, node_count = "S8", 8


class MITC4(Element):
    type_name, node_count = "MITC4", 4


class MITC8(Element):
    type_name, node_count = "MITC8", 8


class QSPT(Element):
    type_name, node_count = "QSPT", 4


class MITC3FRT(Element):
    type_name, node_count = "MITC3FRT", 3


class MITC4FRT(Element):
    type_name, node_count = "MITC4FRT", 4


class MITC6FRT(Element):
    type_name, node_count = "MITC6FRT", 6


class MITC8FRT(Element):
    type_name, node_count = "MITC8FRT", 8


# -----------------------------------------------------------------------------
# One-node point elements
# -----------------------------------------------------------------------------


class MassElement(Element):
    """One-node MASS topology receiving its value from a MassSection."""

    type_name, node_count = "MASS", 1


class RotaryInertiaElement(Element):
    """One-node ROTARYI topology receiving a RotaryInertiaSection."""

    type_name, node_count = "ROTARYI", 1


class SpringElement(Element):
    """One-node SPRING1 topology receiving a SpringSection."""

    type_name, node_count = "SPRING1", 1


ELEMENT_TYPES: dict[str, type[Element]] = {
    item.type_name: item
    for item in (
        C3D4,
        C3D5,
        C3D6,
        C3D8,
        C3D8R,
        C3D10,
        C3D15,
        C3D20,
        C3D20R,
        B33,
        T3,
        T3D2,
        S3,
        S4,
        MITC4,
        S6,
        S8,
        MITC8,
        QSPT,
        MITC3FRT,
        MITC4FRT,
        MITC6FRT,
        MITC8FRT,
        MassElement,
        RotaryInertiaElement,
        SpringElement,
    )
}


class ElementRepository(IdRepository[Element]):
    """Repository of part-local elements keyed by FEMaster element id."""

    def export(self) -> str:
        """Export elements grouped into deterministic TYPE blocks."""

        groups: dict[str, list[Element]] = {}
        for element in self:
            groups.setdefault(element.type_name, []).append(element)

        result: list[str] = []
        for type_name, elements in groups.items():
            result.append(
                block([
                    keyword("ELEMENT", TYPE=type_name),
                    *(element.export() for element in elements),
                ])
            )
        return "\n\n".join(result)


class Surface(NamedObject):
    """Named element-boundary surface definition.

    Each entry references either a local element id or an ELSET name together
    with the FEMaster boundary-side index. The same object is valid in Part and
    Assembly scope; ownership determines where Project exports it.
    """

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.entries: list[tuple[ElementReference, int]] = []

    def add(self, element: ElementReference, side: int) -> "Surface":
        """Append one element or ELSET side reference."""

        self.entries.append((element, int(side)))
        return self

    def export(self) -> str:
        """Export this surface as one SURFACE block."""

        return block([
            keyword("SURFACE", NAME=self.name, TYPE="ELEMENT"),
            *(csv((element, f"S{side}")) for element, side in self.entries),
        ])
