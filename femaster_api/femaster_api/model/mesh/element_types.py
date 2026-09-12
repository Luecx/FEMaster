"""Registry of supported native FEMaster element type names."""

from .elements.c3d4 import C3D4
from .elements.c3d5 import C3D5
from .elements.c3d6 import C3D6
from .elements.c3d8 import C3D8
from .elements.c3d8r import C3D8R
from .elements.c3d10 import C3D10
from .elements.c3d15 import C3D15
from .elements.c3d20 import C3D20
from .elements.c3d20r import C3D20R
from .elements.b33 import B33
from .elements.t3 import T3
from .elements.t3d2 import T3D2
from .elements.s3 import S3
from .elements.s4 import S4
from .elements.s6 import S6
from .elements.s8 import S8
from .elements.mitc4 import MITC4
from .elements.mitc8 import MITC8
from .elements.qspt import QSPT
from .elements.mitc3frt import MITC3FRT
from .elements.mitc4frt import MITC4FRT
from .elements.mitc6frt import MITC6FRT
from .elements.mitc8frt import MITC8FRT
from .elements.mass_element import MassElement
from .elements.rotary_inertia_element import RotaryInertiaElement
from .elements.spring_element import SpringElement


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
