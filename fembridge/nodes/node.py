class Node:
    """
    Represents a geometric node.
    node_id may be None; Nodes.add() will assign the 0-based id automatically.
    """
    def __init__(self, node_id, x: float, y: float, z: float = 0.0):
        self.node_id = node_id
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def to_tuple(self):
        return (self.x, self.y, self.z)

    def to_femaster(self) -> str:
        if self.node_id is None:
            raise RuntimeError("Node has no id assigned.")
        return f"{self.node_id}, {self.x}, {self.y}, {self.z}"

    def to_asami(self) -> str:
        raise NotImplementedError("Node.to_asami not implemented yet.")
