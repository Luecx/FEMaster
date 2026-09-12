"""Constraint-transformation method shared by structural analysis steps.

The values mirror the native FEMaster constraint-method tokens and are kept in a
small utility module because they configure a step rather than representing an
independent model definition.
"""

from __future__ import annotations

from enum import Enum


class ConstraintMethod(Enum):
    """Numerical treatment used for assembled kinematic constraints."""

    NULLSPACE = "NULLSPACE"
    LAGRANGE = "LAGRANGE"
