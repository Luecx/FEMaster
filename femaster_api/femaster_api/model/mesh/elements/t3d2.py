"""Concrete T3D2 element."""

from ..element import Element


class T3D2(Element):
    """Abaqus-compatible two-node truss alias."""

    type_name = "T3D2"
    node_count = 2
