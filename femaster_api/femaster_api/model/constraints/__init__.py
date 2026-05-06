"""Constraint model objects."""

from .connector_constraint import ConnectorConstraint
from .connector_type import ConnectorType
from .constraint_repository import ConstraintRepository
from .coupling_constraint import CouplingConstraint
from .coupling_type import CouplingType
from .rbm_constraint import RBMConstraint
from .tie_constraint import TieConstraint

Constraint = RBMConstraint | CouplingConstraint | ConnectorConstraint | TieConstraint

__all__ = [
    "ConnectorConstraint",
    "ConnectorType",
    "Constraint",
    "ConstraintRepository",
    "CouplingConstraint",
    "CouplingType",
    "RBMConstraint",
    "TieConstraint",
]
