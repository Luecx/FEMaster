"""Base finite-element connectivity object."""

from ..common.format import csv


class Element:
    """Base element preserving its FEMaster ID and local connectivity."""

    type_name = "ELEMENT"
    node_count: int | None = None

    def __init__(self, id: int, nodes: tuple[int, ...] | list[int]) -> None:
        self.id = int(id)
        self.nodes = tuple(int(node) for node in nodes)
        if self.node_count is not None and len(self.nodes) != self.node_count:
            raise ValueError(
                f"{self.type_name} requires {self.node_count} nodes, got {len(self.nodes)}"
            )

    def export(self) -> str:
        return csv((self.id, *self.nodes))
