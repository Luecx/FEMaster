"""Element objects used by the in-memory model."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.nodes import Node

from .element_topology import ElementTopology, FEMASTER_ELEMENT_TYPES


@dataclass(frozen=True, slots=True, eq=False)
class Element:
    """Finite element connectivity storing node objects.

    Elements intentionally do not store FEMaster export IDs. The writer
    resolves element and node IDs through its export context.
    """

    topology: ElementTopology
    nodes: tuple[Node, ...]

    @property
    def femaster_type(self) -> str:
        """Return the FEMaster element type token for this topology."""

        return FEMASTER_ELEMENT_TYPES[self.topology]
