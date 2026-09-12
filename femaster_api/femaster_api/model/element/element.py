"""Base representation of one part-local finite element.

Every concrete FEMaster element type derives from ``Element`` and defines its
native ``TYPE`` token together with the required connectivity size.  The base
class validates that topology immediately, preserves the sparse element ID and
exports only the connectivity row.  Keyword grouping belongs to
``ElementRepository`` because one ``*ELEMENT`` header represents many rows.

The class deliberately does not own sections, materials or regions; those are
separate model concepts linked by semantic names.
"""

from __future__ import annotations

from ..common.format import csv


class Element:
    """Base finite element with persistent ID and fixed node connectivity."""

    type_name = "ELEMENT"
    node_count: int | None = None

    def __init__(self, id: int, nodes: tuple[int, ...] | list[int]) -> None:
        self.id = int(id)
        self.nodes = tuple(int(node) for node in nodes)

        if self.node_count is not None and len(self.nodes) != self.node_count:
            raise ValueError(
                f"{self.type_name} requires {self.node_count} nodes, "
                f"got {len(self.nodes)}"
            )

    def export(self) -> str:
        """Export this element as one connectivity row."""

        return csv((self.id, *self.nodes))
