"""Constraint handling enum."""

from __future__ import annotations

from enum import Enum


class ConstraintMethod(Enum):
    NULLSPACE = "NULLSPACE"
    LAGRANGE = "LAGRANGE"
