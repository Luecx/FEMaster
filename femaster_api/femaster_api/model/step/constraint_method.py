"""Constraint transformation methods."""

from enum import Enum


class ConstraintMethod(Enum):
    """Supported constraint transformation methods."""

    NULLSPACE = "NULLSPACE"
    LAGRANGE = "LAGRANGE"
