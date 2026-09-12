"""Typed named regions used by parts and the compiled assembly.

File names intentionally share the ``region_`` prefix so related concepts remain
adjacent in directory listings: ``region_node``, ``region_element``,
``region_surface`` and ``region_line``.
"""

from .region import Region
from .region_element import ElementRegion
from .region_line import LineRegion
from .region_node import NodeRegion
from .region_repository import RegionRepository
from .region_surface import SurfaceRegion

__all__ = [
    "ElementRegion",
    "LineRegion",
    "NodeRegion",
    "Region",
    "RegionRepository",
    "SurfaceRegion",
]
