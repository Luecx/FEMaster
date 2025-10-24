
from __future__ import annotations
from typing import List, Optional
from .nodeset import NodeSet

class NodeSets:
    def __init__(self) -> None:
        self._items: List[Optional[NodeSet]] = []

    def add(self, s: NodeSet) -> NodeSet:
        self._items.append(s)
        return s

    def to_femaster(self) -> str:
        blocks = []
        for s in self._items:
            if s is None:
                continue
            blocks.append(s.to_femaster())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("NodeSets.to_asami not implemented yet.")
