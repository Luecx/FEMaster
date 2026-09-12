"""Base class for assembly-level non-topological FEMaster features.

Features represent solver definitions that act on compiled entities but are not
themselves mesh topology, sections, loads or constraints.  Every concrete
feature owns its native representation so the repository remains a simple
ordered owner.
"""

from __future__ import annotations


class Feature:
    """Base class for one non-topological model feature."""

    def export(self) -> str:
        raise NotImplementedError
