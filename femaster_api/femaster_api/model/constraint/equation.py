"""Linear multi-point equation."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint
from .equation_term import EquationTerm


class Equation(Constraint):
    """Linear multi-point equation with arbitrary terms."""

    def __init__(self, terms: Iterable[EquationTerm] = ()) -> None:
        self.terms = list(terms)

    def add(self, node: int | str, dof: int, coefficient: float) -> "Equation":
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
