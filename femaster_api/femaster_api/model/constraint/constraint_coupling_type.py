"""Supported coupling formulations for ``Coupling``.

The enum mirrors native FEMaster tokens and keeps formulation selection
type-safe without embedding string literals throughout user code.
"""

from __future__ import annotations

from enum import Enum


class CouplingType(Enum):
    """Kinematic relation used by a coupling constraint."""

    KINEMATIC = "KINEMATIC"
    DISTRIBUTING = "DISTRIBUTING"
