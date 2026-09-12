"""Coupling formulations."""

from enum import Enum


class CouplingType(Enum):
    """Supported coupling formulations."""

    KINEMATIC = "KINEMATIC"
    DISTRIBUTING = "DISTRIBUTING"
