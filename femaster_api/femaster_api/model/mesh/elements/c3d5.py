"""Concrete C3D5 element."""

from ..element import Element


class C3D5(Element):
    """Five-node pyramid accepted by FEMaster and expanded internally to C3D8."""

    type_name = "C3D5"
    node_count = 5
