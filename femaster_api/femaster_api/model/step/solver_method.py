"""Linear solver methods."""

from enum import Enum


class SolverMethod(Enum):
    """Supported linear solver methods."""

    DIRECT = "DIRECT"
    INDIRECT = "INDIRECT"
