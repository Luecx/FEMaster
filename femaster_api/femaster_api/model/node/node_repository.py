"""Repository and grouped export for part-local finite-element nodes.

``NodeRepository`` preserves sparse user IDs through ``IdRepository`` and owns
only the grouping required by the native deck: all node rows of one part are
emitted below a single ``*NODE`` keyword.  It does not renumber, sort or otherwise
reinterpret node identity.
"""

from __future__ import annotations

from ..common.format import block, keyword
from ..common.id_repository import IdRepository
from .node import Node


class NodeRepository(IdRepository[Node]):
    """Own part-local nodes keyed by their FEMaster node ID."""

    def export(self) -> str:
        """Export all owned nodes as one deterministic ``*NODE`` block."""

        if not self._items:
            return ""
        return block([keyword("NODE"), *(node.export() for node in self)])
