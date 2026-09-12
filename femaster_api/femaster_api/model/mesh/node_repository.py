"""Repository for part-local nodes."""

from ..common.format import block, keyword
from ..common.id_repository import IdRepository
from .node import Node


class NodeRepository(IdRepository[Node]):
    """Node repository keyed by semantic FEMaster node ID."""

    def export(self) -> str:
        if not self._items:
            return ""
        return block([keyword("NODE"), *(node.export() for node in self)])
