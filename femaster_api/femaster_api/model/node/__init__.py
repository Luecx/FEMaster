"""Part-local FEMaster node model.

The package separates node storage from element connectivity so each domain can
evolve independently while ``Part`` composes both into a complete reusable mesh.
"""

from .node import Node
from .node_repository import NodeRepository

__all__ = ["Node", "NodeRepository"]
