"""Named FEMaster surface definitions.

Element-side and node-based surfaces are distinct concrete classes.  They can be
owned either by a reusable ``Part`` or by the assembly-level ``Project``.
"""

from .surface import Surface
from .surface_element import ElementSurface
from .surface_node import NodeSurface
from .surface_repository import SurfaceRepository

__all__ = ["ElementSurface", "NodeSurface", "Surface", "SurfaceRepository"]
