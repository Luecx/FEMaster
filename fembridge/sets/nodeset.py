from __future__ import annotations
from typing import Iterable, List, Union, Iterator
from ..nodes.node import Node

NodeListLike = Union[Node, Iterable[Node], "NodeSet"]


class NodeSet:
    """
    Container nur für Node-Objekte (niemals ints).
    name kann None sein (interne/temporäre Sets).
    """

    def __init__(self, name: str | None = None, nodes: List[Node] | None = None) -> None:
        self.name: str | None = (str(name) if name else None)
        self.nodes: List[Node] = []
        if nodes:
            for n in nodes:
                if not isinstance(n, Node):
                    raise TypeError("NodeSet darf nur Node-Instanzen enthalten.")
            self.nodes.extend(nodes)

    # ---------------- Fabriken ----------------

    @staticmethod
    def internal(value: NodeListLike) -> NodeSet:
        """Unbenanntes Set aus Node | Iterable[Node] | NodeSet."""
        if isinstance(value, NodeSet):
            return value
        return NodeSet(None, NodeSet._normalize(value))

    @staticmethod
    def create(name: str, value: NodeListLike) -> NodeSet:
        """Benanntes Set (kopiert Nodes falls NodeSet)."""
        if not name:
            raise ValueError("NodeSet.create(): name darf nicht leer sein.")
        if isinstance(value, NodeSet):
            return NodeSet(name, list(value.nodes))
        return NodeSet(name, NodeSet._normalize(value))

    # --------------- Normalisierung ---------------

    @staticmethod
    def _normalize(value: NodeListLike) -> List[Node]:
        """Konvertiert Node | Iterable[Node] | NodeSet -> list[Node]."""
        if isinstance(value, Node):
            return [value]
        if isinstance(value, NodeSet):
            return list(value.nodes)
        try:
            lst = list(value)  # type: ignore[arg-type]
        except TypeError:
            raise TypeError("Erwarte Node, Iterable[Node] oder NodeSet.")
        for n in lst:
            if not isinstance(n, Node):
                raise TypeError("Iterable muss ausschließlich Node-Instanzen enthalten.")
        return lst

    # ---------------- Convenience ----------------

    def add(self, value: NodeListLike) -> None:
        """Fügt Node | Iterable[Node] | NodeSet hinzu."""
        self.nodes.extend(NodeSet._normalize(value))

    def __iter__(self) -> Iterator[Node]:
        return iter(self.nodes)

    def __len__(self) -> int:
        return len(self.nodes)

    def __bool__(self) -> bool:
        return bool(self.nodes)

    # ---------------- Mengen-Operatoren ----------------

    def __or__(self, other: NodeSet) -> NodeSet:
        """Vereinigung (stabile Reihenfolge, Deduplizierung nach node.node_id)."""
        if not isinstance(other, NodeSet):
            return NotImplemented
        out: List[Node] = []
        seen: set[int] = set()
        for n in self.nodes + other.nodes:
            nid = n.node_id
            if nid not in seen:
                out.append(n)
                seen.add(nid)
        return NodeSet(None, out)

    def __and__(self, other: NodeSet) -> NodeSet:
        """Schnittmenge (Reihenfolge nach self, Deduplizierung nach node.node_id)."""
        if not isinstance(other, NodeSet):
            return NotImplemented
        other_ids = {n.node_id for n in other.nodes}
        out: List[Node] = []
        seen: set[int] = set()
        for n in self.nodes:
            nid = n.node_id
            if nid in other_ids and nid not in seen:
                out.append(n)
                seen.add(nid)
        return NodeSet(None, out)

    def __sub__(self, other: "NodeSet") -> "NodeSet":
        """Mengendifferenz: alle Nodes aus self, die NICHT in other sind."""
        if not isinstance(other, NodeSet):
            return NotImplemented
        other_ids = {n.node_id for n in other.nodes}
        out: List[Node] = []
        seen: set[int] = set()
        for n in self.nodes:
            nid = n.node_id
            if nid not in other_ids and nid not in seen:
                out.append(n)
                seen.add(nid)
        return NodeSet(None, out)

    # ---------------- Exports ----------------

    def to_femaster(self) -> str:
        if not self.name:
            raise ValueError("Unbenanntes NodeSet kann nicht exportiert werden. Verwende NodeSet.create(name, ...).")

        ids = [n.node_id for n in self.nodes]
        lines = [f"*NSET, NAME={self.name}"]
        line: List[str] = []
        for i, nid in enumerate(ids, 1):
            line.append(str(nid))
            if i % 16 == 0:
                lines.append(", ".join(line))
                line = []
        if line:
            lines.append(", ".join(line))
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError
