"""Solver execution devices."""

from enum import Enum


class SolverDevice(Enum):
    """Supported solver execution devices."""

    CPU = "CPU"
    GPU = "GPU"
