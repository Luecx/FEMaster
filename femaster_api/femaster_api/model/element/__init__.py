"""Part-local finite-element topology.

Every supported native element formulation has its own module and concrete
class.  ``ElementRepository`` owns sparse IDs and deck grouping, while the
individual classes carry only formulation-specific type and connectivity data.
"""

from .element import Element
from .element_repository import ElementRepository
from .element_types import ELEMENT_TYPES
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

__all__ = ['Element', 'ElementRepository', 'ELEMENT_TYPES', 'C3D4', 'C3D5', 'C3D6', 'C3D8', 'C3D8R', 'C3D10', 'C3D15', 'C3D20', 'C3D20R', 'B33', 'T3', 'T3D2', 'S3', 'S4', 'S6', 'S8', 'MITC4', 'MITC8', 'QSPT', 'MITC3FRT', 'MITC4FRT', 'MITC6FRT', 'MITC8FRT', 'MassElement', 'RotaryInertiaElement', 'SpringElement']
