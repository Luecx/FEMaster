"""One coefficient term of a linear multi-point equation.

An equation term identifies one node reference and one structural degree of
freedom together with its scalar coefficient.  The small value object is kept in
its own module to preserve the project's strict one-class-per-file rule.
"""

from __future__ import annotations


class EquationTerm:
    """One ``coefficient * u(node, dof)`` contribution."""

    def __init__(self, node: int | str, dof: int, coefficient: float) -> None:
        self.node = node
        self.dof = int(dof)
        self.coefficient = float(coefficient)
        if self.dof < 1 or self.dof > 6:
            raise ValueError("EquationTerm dof must be between 1 and 6")
