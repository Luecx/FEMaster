from typing import Sequence, Union
from ..nodes.node import Node

NodeLike = Union[int, Node]

class Element:
    """
    Base finite element class.

    Accepts either node IDs or Node objects in the constructor.
    If Node objects are provided, they must already have a node_id assigned
    (via Nodes.add or manually). Otherwise, an error is raised.

    Example:
        Element(0, [n1, n2], "B33")
        Element(1, [0, 1], "B33")
    """
    def __init__(self, element_id, node_refs: Sequence[NodeLike], elem_type: str | None = None):
        def _to_id(node: NodeLike) -> int:
            if isinstance(node, Node):
                if node.node_id is None:
                    raise ValueError("Node object must have a node_id before being used in an Element.")
                return int(node.node_id)
            return int(node)

        self.element_id = element_id
        self.node_ids = [_to_id(n) for n in node_refs]
        self.elem_type = elem_type or self.__class__.__name__

    # --- geometric helpers ---
    def to_second_order(self, new_node_ids):
        raise NotImplementedError("Should be implemented by subclass")

    def compute_center(self, geometry):
        coords = [geometry.nodes[nid].to_tuple() for nid in self.node_ids]
        x = sum(c[0] for c in coords) / len(coords)
        y = sum(c[1] for c in coords) / len(coords)
        z = sum(c[2] for c in coords) / len(coords)
        return (x, y, z)

    def connectivity(self):
        raise NotImplementedError("Should be implemented by subclass")

    def subdivide(self, edge_nodes, geometry, only_quads=False):
        raise NotImplementedError("Should be implemented by subclass")

    def mirror_ids(self):
        raise NotImplementedError("Should be implemented by subclass")

    # --- serialization ---
    def femaster_record(self) -> str:
        """Return the line 'id, n1, n2, ...' for FEMaster export."""
        if self.element_id is None:
            raise RuntimeError("Element has no id assigned.")
        return f"{self.element_id}, " + ", ".join(str(i) for i in self.node_ids)

    def to_asami(self) -> str:
        raise NotImplementedError(f"to_asami not implemented yet for {self.elem_type}.")
