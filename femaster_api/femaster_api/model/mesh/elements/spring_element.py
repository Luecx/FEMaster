"""Concrete SPRING1 element."""

from ..element import Element


class SpringElement(Element):
    """One-node SPRING1 topology receiving a SpringSection."""

    type_name = "SPRING1"
    node_count = 1
