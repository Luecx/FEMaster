"""Concrete T3 element."""

from ..element import Element


class T3(Element):
    """Two-node native truss type."""

    type_name = "T3"
    node_count = 2
