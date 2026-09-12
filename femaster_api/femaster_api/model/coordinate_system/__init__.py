"""Global coordinate-system definitions used across the model.

Rectangular and cylindrical coordinate systems are distinct concrete classes
because their defining geometry and native syntax differ.  The shared named
repository allows loads, supports and sections to reference either form through
stable semantic names.
"""

from .coordinate_system import CoordinateSystem
from .coordinate_system_cylindrical import CylindricalCoordinateSystem
from .coordinate_system_rectangular import RectangularCoordinateSystem
from .coordinate_system_repository import CoordinateSystemRepository

__all__ = [
    "CoordinateSystem",
    "CoordinateSystemRepository",
    "CylindricalCoordinateSystem",
    "RectangularCoordinateSystem",
]
