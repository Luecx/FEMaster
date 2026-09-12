"""Solver execution device used by analysis steps.

The enum mirrors FEMaster's native ``DEVICE`` tokens.  It lives in ``step/util``
because device selection is an analysis control shared by several concrete step
types rather than a top-level model entity.
"""

from __future__ import annotations

from enum import Enum


class SolverDevice(Enum):
    """Execution device selected for one linear-system solve."""

    CPU = "CPU"
    GPU = "GPU"
