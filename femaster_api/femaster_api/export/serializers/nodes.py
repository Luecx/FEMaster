"""FEMaster serialization for nodes."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.nodes import NodeRepository


def write_nodes(nodes: NodeRepository) -> str:
    if len(nodes) == 0:
        return ""
    lines = [keyword("NODE")]
    lines.extend(csv((node.id, node.x, node.y, node.z)) for node in nodes.all())
    return block(lines)
