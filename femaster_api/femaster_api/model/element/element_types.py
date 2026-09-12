"""Registry of canonical concrete element classes understood by the INP reader.

Only FEMaster's primary internal element names are registered here.  Historical
aliases and implementation-specific shell/point variants are deliberately not
promoted to first-class Python element types.  Unknown ``TYPE`` tokens can still
be preserved by ``Project.read_inp`` through the generic ``Element`` fallback,
so reducing this mapping does not make the input reader destructive.
"""

from __future__ import annotations

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


ELEMENT_TYPES = {
    "B33": B33,
    "C3D4": C3D4,
    "C3D5": C3D5,
    "C3D6": C3D6,
    "C3D8": C3D8,
    "C3D8R": C3D8R,
    "C3D10": C3D10,
    "C3D15": C3D15,
    "C3D20": C3D20,
    "C3D20R": C3D20R,
    "S3": S3,
    "S4": S4,
    "S6": S6,
    "S8": S8,
    "T3": T3,
}
