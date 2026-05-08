"""Independent VTK mesh export helpers for FEMaster models."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model import Element, ElementTopology, Model, Node


VTK_CELL_TYPES: dict[ElementTopology, int] = {
    ElementTopology.TRUSS2: 3,
    ElementTopology.BEAM2: 3,
    ElementTopology.TRI3: 5,
    ElementTopology.QUAD4: 9,
    ElementTopology.MITC4: 9,
    ElementTopology.QSPT: 9,
    ElementTopology.TET4: 10,
    ElementTopology.HEX8: 12,
    ElementTopology.WEDGE6: 13,
    ElementTopology.PYRAMID5: 14,
    ElementTopology.TRI6: 22,
    ElementTopology.QUAD8: 23,
    ElementTopology.TET10: 24,
    ElementTopology.HEX20: 25,
    ElementTopology.HEX20R: 25,
    ElementTopology.WEDGE15: 26,
}


@dataclass(frozen=True, slots=True)
class VtkMesh:
    """Flat unstructured-grid arrays derived from a model."""

    points: tuple[tuple[float, float, float], ...]
    connectivity: tuple[int, ...]
    offsets: tuple[int, ...]
    cell_types: tuple[int, ...]
    node_ids: dict[Node, int]
    element_ids: dict[Element, int]


def build_vtk_mesh(model: Model) -> VtkMesh:
    """Create VTK unstructured-grid arrays using zero-based node indices."""

    node_ids = {node: index for index, node in enumerate(model.nodes)}
    element_ids = {element: index for index, element in enumerate(model.elements)}

    points: list[tuple[float, float, float]] = []
    for node in model.nodes:
        points.append((node.x, node.y, node.z))

    connectivity: list[int] = []
    offsets: list[int] = []
    cell_types: list[int] = []

    for element in model.elements:
        try:
            cell_type = VTK_CELL_TYPES[element.topology]
        except KeyError as exc:
            raise ValueError(f"unsupported VTK topology: {element.topology.value}") from exc
        connectivity.extend(node_ids[node] for node in element.nodes)
        offsets.append(len(connectivity))
        cell_types.append(cell_type)

    return VtkMesh(
        points=tuple(points),
        connectivity=tuple(connectivity),
        offsets=tuple(offsets),
        cell_types=tuple(cell_types),
        node_ids=node_ids,
        element_ids=element_ids,
    )
