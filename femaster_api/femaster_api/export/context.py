"""Export-only ID maps for FEMaster deck serialization."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.elements import Element
from femaster_api.model.model import Model
from femaster_api.model.nodes import Node
from femaster_api.model.surfaces import Surface


@dataclass(frozen=True)
class ExportContext:
    """Local FEMaster ID maps built only while writing a deck."""

    node_ids: dict[Node, int]
    element_ids: dict[Element, int]
    surface_ids: dict[Surface, int]

    @classmethod
    def from_model(cls, model: Model) -> "ExportContext":
        """Build one-based export IDs from repository order."""

        return cls(
            node_ids={node: index for index, node in enumerate(model.nodes, start=1)},
            element_ids={element: index for index, element in enumerate(model.elements, start=1)},
            surface_ids={surface: index for index, surface in enumerate(model.surfaces, start=1)},
        )

    def node_id(self, node: Node) -> int:
        return self.node_ids[node]

    def element_id(self, element: Element) -> int:
        return self.element_ids[element]

    def surface_id(self, surface: Surface) -> int:
        return self.surface_ids[surface]

    def target_token(self, target: object) -> str:
        """Return the FEMaster token for a model object or named target."""

        name = getattr(target, "name", None)
        if isinstance(name, str) and name:
            return name
        if isinstance(target, Node):
            return str(self.node_id(target))
        if isinstance(target, Element):
            return str(self.element_id(target))
        if isinstance(target, Surface):
            return str(self.surface_id(target))
        return str(target)
