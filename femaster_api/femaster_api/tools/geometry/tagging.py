"""Mesh tagging helpers."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from femaster_api.model.elements import Element
from femaster_api.model.model import Model
from femaster_api.model.nodes import Node
from femaster_api.model.sets import ElementSet, NodeSet

from femaster_api.tools.geometry.base import GeometryEntity
from .vector import Point3, centroid


def apply_geometry_tags(
    model: Model,
    entities: Iterable[GeometryEntity],
    *,
    tolerance: float = 1.0e-6,
) -> None:
    """Create node and element sets from geometry tags.

    Nodes are classified by their coordinates. Elements are classified by the
    centroid of their connected nodes.
    """

    node_members: dict[str, list[Node]] = defaultdict(list)
    element_members: dict[str, list[Element]] = defaultdict(list)
    tagged_entities = tuple(entity for entity in entities if entity.tags)

    for node in model.nodes:
        point = _node_point(node)
        for entity in tagged_entities:
            if entity.contains(point, tolerance=tolerance):
                for tag in entity.tags:
                    node_members[tag].append(node)

    for element in model.elements:
        center = centroid(_node_point(node) for node in element.nodes)
        for entity in tagged_entities:
            if entity.contains(center, tolerance=tolerance):
                for tag in entity.tags:
                    element_members[tag].append(element)

    for tag, nodes in node_members.items():
        if nodes:
            model.sets.add(NodeSet(tag, tuple(dict.fromkeys(nodes))))
    for tag, elements in element_members.items():
        if elements:
            model.sets.add(ElementSet(tag, tuple(dict.fromkeys(elements))))


def _node_point(node: Node) -> Point3:
    return (float(node.x), float(node.y), float(node.z))
