from __future__ import annotations
from typing import List, Optional
from .node import Node

class Nodes:
    """
    Sparse 0-based list of Node objects.
    IDs are assigned automatically if Node.node_id is None.
    """
    def __init__(self) -> None:
        self._items: List[Optional[Node]] = []

    def add(self, node: Node) -> Node:
        if node.node_id is None:
            for i, v in enumerate(self._items):
                if v is None:
                    node.node_id = i
                    self._items[i] = node
                    return node
            node.node_id = len(self._items)
            self._items.append(node)
            return node
        else:
            nid = int(node.node_id)
            if nid < 0:
                raise ValueError("node_id must be >= 0")
            if nid >= len(self._items):
                self._items.extend([None] * (nid + 1 - len(self._items)))
            self._items[nid] = node
            return node

    def __getitem__(self, node_id: int) -> Node:
        val = self._items[node_id]
        if val is None:
            raise KeyError(f"Node id {node_id} is empty (deleted or uninitialized).")
        return val

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        lines = ["*NODE"]
        for n in self._items:
            if n is None:
                continue
            lines.append(n.to_femaster())
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError("Nodes.to_asami not implemented yet.")
