"""FEMaster serialization for elements."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.elements import ElementRepository, ElementTopology


def write_elements(elements: ElementRepository, context: ExportContext) -> str:
    if len(elements) == 0:
        return ""

    grouped: dict[ElementTopology, list[object]] = {}
    for element in elements:
        grouped.setdefault(element.topology, []).append(element)

    lines: list[str] = []
    for topology, group in sorted(grouped.items(), key=lambda item: item[0].value):
        if lines:
            lines.append("")
        sample = group[0]
        lines.append(keyword("ELEMENT", TYPE=sample.femaster_type))
        for element in group:
            node_ids = tuple(context.node_id(node) for node in element.nodes)
            lines.append(csv((context.element_id(element), *node_ids)))
    return block(lines)
