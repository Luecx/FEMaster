"""Linear-system solution strategy used by analysis steps.

``SolverMethod`` is a shared numerical control and therefore belongs in
``step/util`` rather than beside concrete loadcase classes.
"""

from __future__ import annotations

from enum import Enum


class SolverMethod(Enum):
    """High-level direct or iterative solver selection."""

    DIRECT = "DIRECT"
    INDIRECT = "INDIRECT"
