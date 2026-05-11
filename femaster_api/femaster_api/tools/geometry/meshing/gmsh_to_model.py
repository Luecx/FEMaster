"""Convert Gmsh meshes to model repositories."""

from __future__ import annotations

from femaster_api.model.elements import Element, ElementTopology
from femaster_api.model.model import Model
from femaster_api.model.nodes import Node


def append_gmsh_mesh(model: Model, gmsh, dimension: int) -> None:
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    node_by_tag: dict[int, Node] = {}
    for index, tag in enumerate(node_tags):
        offset = index * 3
        node_by_tag[int(tag)] = model.nodes.add(Node(float(coords[offset]), float(coords[offset + 1]), float(coords[offset + 2])))

    element_types, _, element_node_tags = gmsh.model.mesh.getElements(dimension)
    for element_type, connectivity in zip(element_types, element_node_tags):
        topology = _topology_for_gmsh_type(int(element_type))
        width = _gmsh_node_count(int(element_type))
        for start in range(0, len(connectivity), width):
            nodes = tuple(node_by_tag[int(tag)] for tag in connectivity[start : start + width])
            model.elements.add(Element(topology, nodes))


def _topology_for_gmsh_type(element_type: int) -> ElementTopology:
    mapping = {
        1: ElementTopology.BEAM2,
        2: ElementTopology.TRI3,
        3: ElementTopology.QUAD4,
        4: ElementTopology.TET4,
        5: ElementTopology.HEX8,
        6: ElementTopology.WEDGE6,
        7: ElementTopology.PYRAMID5,
        8: ElementTopology.BEAM2,
        9: ElementTopology.TRI6,
        11: ElementTopology.TET10,
        16: ElementTopology.QUAD8,
    }
    try:
        return mapping[element_type]
    except KeyError as exc:
        raise ValueError(f"unsupported gmsh element type: {element_type}") from exc


def _gmsh_node_count(element_type: int) -> int:
    counts = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 11: 10, 16: 8}
    try:
        return counts[element_type]
    except KeyError as exc:
        raise ValueError(f"unsupported gmsh element type: {element_type}") from exc
