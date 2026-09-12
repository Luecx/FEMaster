"""Assembly-level FEMaster constraint definitions.

Concrete constraint files use the ``constraint_*`` prefix and own their native
keyword syntax.  ``ConstraintRepository`` is intentionally heterogeneous and
ordered: it preserves model definition order without introducing a central
serializer switch over constraint types.

Small helper/value types that are meaningful only inside one constraint remain
nested on their owning class.  In particular, connector formulations are exposed
as ``Connector.Type`` and equation terms as ``Equation.Term`` rather than as
separate package-level model concepts.
"""

from .constraint import Constraint
from .constraint_connector import Connector
from .constraint_coupling import Coupling
from .constraint_coupling_type import CouplingType
from .constraint_equation import Equation
from .constraint_repository import ConstraintRepository
from .constraint_rigid_body import RigidBodyConstraint
from .constraint_tie import Tie

__all__ = [
    "Connector",
    "Constraint",
    "ConstraintRepository",
    "Coupling",
    "CouplingType",
    "Equation",
    "RigidBodyConstraint",
    "Tie",
]
