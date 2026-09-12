"""Part-local finite-element node."""

from ..common.format import csv


class Node:
    """Part-local node preserving its FEMaster ID."""

    def __init__(self, id: int, x: float, y: float, z: float = 0.0) -> None:
        self.id = int(id)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def export(self) -> str:
        return csv((self.id, self.x, self.y, self.z))
