"""Linear multi-point equation whose terms reference real node-domain objects.

``Equation`` owns the complete native ``*EQUATION`` relation.  Each nested
``Equation.Term`` stores a ``Node`` or ``NodeRegion`` object, never an integer ID
or semantic-name string.  The object is converted to its native token only when
serializing the term triple.

``Term`` remains nested because it has no independent model lifecycle.  Keeping
it named still gives callers readable attributes while avoiding a second public
constraint concept.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..node.node import Node
from ..region.region_node import NodeRegion
from .constraint import Constraint


class Equation(Constraint):
    """Linear relation between nodal degrees of freedom."""

    class Term:
        """One ``coefficient * u(node, dof)`` contribution to an equation."""

        def __init__(
            self,
            node: Node | NodeRegion,
            dof: int,
            coefficient: float,
        ) -> None:
            if not isinstance(node, (Node, NodeRegion)):
                raise TypeError("node must be Node or NodeRegion")
            self.node = node
            self.dof = int(dof)
            self.coefficient = float(coefficient)

            if self.dof < 1 or self.dof > 6:
                raise ValueError("Equation term dof must be between 1 and 6")

    def __init__(
        self,
        terms: Iterable[Term | tuple[Node | NodeRegion, int, float]] = (),
    ) -> None:
        self.terms: list[Equation.Term] = []

        for term in terms:
            if isinstance(term, Equation.Term):
                self.add(term.node, term.dof, term.coefficient)
            else:
                node, dof, coefficient = term
                self.add(node, dof, coefficient)

    def add(
        self,
        node: Node | NodeRegion,
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
            node = term.node.id if isinstance(term.node, Node) else term.node.name
            values.extend((node, term.dof, term.coefficient))

        return block([
            keyword("EQUATION"),
            csv((len(self.terms),)),
            csv(values),
        ])
