"""Node-region specialization storing actual ``Node`` instances.

``NodeRegion`` is the object-model counterpart of native ``*NSET``.  Membership
is strongly typed as ``Node`` objects; sparse node IDs are only recovered at
export time.  This lets loads, supports and constraints reference the region
object directly without introducing a second string/ID reference layer.
"""

from __future__ import annotations

from ..node.node import Node
from .region import Region


class NodeRegion(Region[Node]):
    """Ordered region of concrete nodes exported through ``*NSET``."""

    keyword_name = "NSET"
    name_key = "NAME"
    member_type = Node

    def _member_value(self, member: Node) -> int:
        return member.id
