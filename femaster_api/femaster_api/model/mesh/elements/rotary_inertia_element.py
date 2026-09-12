"""Concrete ROTARYI element."""

from ..element import Element


class RotaryInertiaElement(Element):
    """One-node ROTARYI topology receiving a RotaryInertiaSection."""

    type_name = "ROTARYI"
    node_count = 1
