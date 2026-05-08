"""FEMaster serialization for nodes."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.nodes import NodeRepository


def write_nodes(nodes: NodeRepository, context: ExportContext) -> str:
    if len(nodes) == 0:
        return ""
    lines = [keyword("NODE")]
    lines.extend(csv((context.node_id(node), node.x, node.y, node.z)) for node in nodes)
    return block(lines)
