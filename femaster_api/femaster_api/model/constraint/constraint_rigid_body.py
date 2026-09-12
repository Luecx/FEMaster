"""Rigid-body-mode constraint/removal on an element region.

The object maps directly to FEMaster's ``*RBM`` keyword and stores only the
semantic element-region reference required by that operation.
"""

from __future__ import annotations

from ..common.format import keyword
from .constraint import Constraint


class RigidBodyConstraint(Constraint):
    """Rigid-body constraint acting on one element region."""

    def __init__(self, element_region: str) -> None:
        self.element_region = str(element_region)

    def export(self) -> str:
        return keyword("RBM", ELSET=self.element_region)
