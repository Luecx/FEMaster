"""Base class for assembly-level FEMaster constraints.

Constraints connect or restrict already-defined assembly entities.  Every
concrete constraint owns its native keyword representation so the repository can
remain a simple ordered heterogeneous container without type-switching
serialization logic.
"""

from __future__ import annotations


class Constraint:
    """Base class for one assembly-level kinematic constraint."""

    def export(self) -> str:
        raise NotImplementedError
