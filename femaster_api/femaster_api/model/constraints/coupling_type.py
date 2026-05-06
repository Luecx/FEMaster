"""Coupling constraint type enum."""

from __future__ import annotations

from enum import Enum


class CouplingType(Enum):
    KINEMATIC = "KINEMATIC"
    STRUCTURAL = "STRUCTURAL"
