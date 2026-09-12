"""Linear multi-point equation with an arbitrary number of terms.

The class owns equation-term order because FEMaster's input form begins with the
term count followed by node/DOF/coefficient triples.  Terms can be appended
fluently without introducing a separate builder object.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint
from .constraint_equation_term import EquationTerm


class Equation(Constraint):
    """Linear relation between nodal degrees of freedom."""

    def __init__(self, terms: Iterable[EquationTerm] = ()) -> None:
        self.terms = list(terms)

    def add(
        self,
        node: int | str,
        dof: int,
        coefficient: float,
    ) -> "Equation":
        """Append one equation term and return this equation."""

        self.terms.append(EquationTerm(node, dof, coefficient))
        return self

    def export(self) -> str:
        values: list[object] = []
        for term in self.terms:
            values.extend((term.node, term.dof, term.coefficient))

        return block([
            keyword("EQUATION"),
            csv((len(self.terms),)),
            csv(values),
        ])
