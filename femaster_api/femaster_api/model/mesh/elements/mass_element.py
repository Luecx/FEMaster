"""Concrete MASS element."""

from ..element import Element


class MassElement(Element):
    """One-node MASS topology receiving a MassSection."""

    type_name = "MASS"
    node_count = 1
