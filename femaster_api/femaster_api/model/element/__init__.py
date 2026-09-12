"""Canonical finite-element classes exposed by the Python model.

The Python API intentionally mirrors FEMaster's primary internal element names
instead of accumulating aliases or implementation-specific formulations.  Each
concrete element lives in its own module; this package only re-exports the
supported canonical solids, shells, truss and beam together with the common
``Element`` base class and sparse-ID repository.
"""

from .element import Element
from .element_repository import ElementRepository
from .element_types import ELEMENT_TYPES
from .element_b33 import B33
from .element_c3d4 import C3D4
from .element_c3d5 import C3D5
from .element_c3d6 import C3D6
from .element_c3d8 import C3D8
from .element_c3d8r import C3D8R
from .element_c3d10 import C3D10
from .element_c3d15 import C3D15
from .element_c3d20 import C3D20
from .element_c3d20r import C3D20R
from .element_s3 import S3
from .element_s4 import S4
from .element_s6 import S6
from .element_s8 import S8
from .element_t3 import T3

__all__ = [
    "Element",
    "ElementRepository",
    "ELEMENT_TYPES",
    "B33",
    "C3D4",
    "C3D5",
    "C3D6",
    "C3D8",
    "C3D8R",
    "C3D10",
    "C3D15",
    "C3D20",
    "C3D20R",
    "S3",
    "S4",
    "S6",
    "S8",
    "T3",
]
