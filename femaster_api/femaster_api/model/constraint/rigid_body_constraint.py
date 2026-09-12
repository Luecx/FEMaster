"""Rigid-body constraint."""

from ..common.format import keyword
from .constraint import Constraint


class RigidBodyConstraint(Constraint):
    """Rigid-body removal/constraint acting on an element region."""

    def __init__(self, element_region: str) -> None:
        self.element_region = str(element_region)

    def export(self) -> str:
        return keyword("RBM", ELSET=self.element_region)
