"""Central model container for the independent FEMaster Python API."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from .analysis import StepRepository
from .constraints import ConstraintRepository
from .elements import Element, ElementRepository, ElementTopology
from .features import FeatureRepository
from .fields import FieldRepository
from .loads import LoadCollectorRepository, LoadRepository
from .materials import MaterialRepository
from .nodes import Node, NodeRepository
from .orientations import OrientationRepository
from .sections import SectionRepository
from .sets import ElementSet, NodeSet, SetRepository
from .supports import SupportCollectorRepository, SupportRepository
from .surfaces import SurfaceRepository


class Model:
    """Neutral finite-element model container.

    The model does not render FEMaster text itself. Use
    `femaster_api.export.FEMasterWriter` for deck generation.
    """

    def __init__(self, name: str = "model") -> None:
        self.name = name
        self.nodes = NodeRepository()
        self.elements = ElementRepository()
        self.sets = SetRepository()
        self.surfaces = SurfaceRepository()
        self.orientations = OrientationRepository()
        self.materials = MaterialRepository()
        self.sections = SectionRepository()
        self.fields = FieldRepository()
        self.features = FeatureRepository()
        self.loads = LoadRepository()
        self.load_collectors = LoadCollectorRepository()
        self.supports = SupportRepository()
        self.support_collectors = SupportCollectorRepository()
        self.constraints = ConstraintRepository()
        self.steps = StepRepository()

    def add_node(self, x: float, y: float, z: float) -> Node:
        return self.nodes.add(Node(x, y, z))

    def add_element(
        self,
        topology: ElementTopology,
        *nodes: Node,
    ) -> Element:
        connectivity = _normalize_connectivity(nodes)
        return self.elements.add(Element(topology, connectivity))

    def add_element_from_nodes(
        self,
        topology: ElementTopology,
        nodes: Iterable[Node],
        **kwargs: object,
    ) -> Element:
        return self.add_element(topology, *tuple(nodes), **kwargs)

    def add_material(self, material):
        return self.materials.add(material)

    def add_section(self, section):
        return self.sections.add(section)

    def create_node_set(self, name: str, nodes: Iterable[object]):
        return self.sets.add(NodeSet(name, tuple(nodes)))

    def create_element_set(self, name: str, elements: Iterable[object]):
        return self.sets.add(ElementSet(name, tuple(elements)))

    def add_step(self, step):
        return self.steps.add(step)

    def validate(self, *, raise_on_error: bool = False):
        from femaster_api.validation.checks import validate_model

        diagnostics = validate_model(self)
        if raise_on_error:
            diagnostics.raise_for_errors()
        return diagnostics

    def write_deck(self, path: str | Path) -> None:
        """Compatibility convenience around FEMasterWriter."""

        from femaster_api.export.femaster_writer import FEMasterWriter

        FEMasterWriter(self).write(path)


def _normalize_connectivity(nodes: tuple[Node, ...]) -> tuple[Node, ...]:
    if len(nodes) == 1 and not isinstance(nodes[0], Node):
        value = nodes[0]
        try:
            return tuple(value)  # type: ignore[arg-type]
        except TypeError:
            pass
    return nodes
