"""Registry of concrete element classes understood by the INP reader.

The mapping is data, not a serializer dispatch mechanism.  It exists solely so
``Project.read_inp`` can construct the concrete Python class named by an
``*ELEMENT, TYPE=...`` block.  Export remains polymorphic through each element's
own ``type_name`` and ``export`` behavior.
"""

from __future__ import annotations

from .element_c3d4 import C3D4
from .element_c3d5 import C3D5
from .element_c3d6 import C3D6
from .element_c3d8 import C3D8
from .element_c3d8r import C3D8R
from .element_c3d10 import C3D10
from .element_c3d15 import C3D15
from .element_c3d20 import C3D20
from .element_c3d20r import C3D20R
from .element_b33 import B33
from .element_t3 import T3
from .element_t3d2 import T3D2
from .element_s3 import S3
from .element_s4 import S4
from .element_s6 import S6
from .element_s8 import S8
from .element_mitc4 import MITC4
from .element_mitc8 import MITC8
from .element_qspt import QSPT
from .element_mitc3frt import MITC3FRT
from .element_mitc4frt import MITC4FRT
from .element_mitc6frt import MITC6FRT
from .element_mitc8frt import MITC8FRT
from .element_mass import MassElement
from .element_rotary_inertia import RotaryInertiaElement
from .element_spring import SpringElement

ELEMENT_TYPES = {
    "C3D4": C3D4,
    "C3D5": C3D5,
    "C3D6": C3D6,
    "C3D8": C3D8,
    "C3D8R": C3D8R,
    "C3D10": C3D10,
    "C3D15": C3D15,
    "C3D20": C3D20,
    "C3D20R": C3D20R,
    "B33": B33,
    "T3": T3,
    "T3D2": T3D2,
    "S3": S3,
    "S4": S4,
    "S6": S6,
    "S8": S8,
    "MITC4": MITC4,
    "MITC8": MITC8,
    "QSPT": QSPT,
    "MITC3FRT": MITC3FRT,
    "MITC4FRT": MITC4FRT,
    "MITC6FRT": MITC6FRT,
    "MITC8FRT": MITC8FRT,
    "MASS": MassElement,
    "ROTARYI": RotaryInertiaElement,
    "SPRING1": SpringElement,
}
