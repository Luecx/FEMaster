"""Linear multi-point equation with an arbitrary number of terms.

``Equation`` owns the complete representation of one native ``*EQUATION``
constraint, including the lightweight value object used for individual terms.
A term has no independent lifecycle in the FEMaster model, so it is deliberately
nested as ``Equation.Term`` instead of being exported as a separate top-level
model class.  This keeps the public constraint namespace compact while retaining
named attributes for node reference, degree of freedom and coefficient.

All terms, including tuples supplied to the constructor, pass through ``add``.
That gives constructor input and fluent builder input exactly the same validation
and normalization rules while preserving the user-defined term order required by
FEMaster's term-count-plus-triples input syntax.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint


class Equation(Constraint):
    """Linear relation between nodal degrees of freedom."""

    class Term:
        """One ``coefficient * u(node, dof)`` contribution to an equation."""

        def __init__(
            self,
            node: int | str,
            dof: int,
            coefficient: float,
        ) -> None:
            self.node = node
            self.dof = int(dof)
            self.coefficient = float(coefficient)

            if self.dof < 1 or self.dof > 6:
                raise ValueError("Equation term dof must be between 1 and 6")

    def __init__(
        self,
        terms: Iterable[Term | tuple[int | str, int, float]] = (),
    ) -> None:
        self.terms: list[Equation.Term] = []

        # Normalize every constructor term through the same path used by add().
        # Existing Term objects are copied so an Equation always owns its term
        # instances and later external mutation cannot change the equation.
        for term in terms:
            if isinstance(term, Equation.Term):
                self.add(term.node, term.dof, term.coefficient)
            else:
                node, dof, coefficient = term
                self.add(node, dof, coefficient)

    def add(
        self,
        node: int | str,
        dof: int,
        coefficient: float,
    ) -> "Equation":
        """Append one validated term and return this equation for fluent use."""

        self.terms.append(self.Term(node, dof, coefficient))
        return self

    def export(self) -> str:
        """Export FEMaster's term count followed by node/DOF/coefficient triples."""

        values: list[object] = []
        for term in self.terms:
            values.extend((term.node, term.dof, term.coefficient))

        return block([
            keyword("EQUATION"),
            csv((len(self.terms),)),
            csv(values),
        ])
