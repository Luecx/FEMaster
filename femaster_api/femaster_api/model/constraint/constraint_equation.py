"""Linear multi-point equation with an arbitrary number of terms.

``Equation`` owns the complete representation of one native ``*EQUATION``
constraint.  A term has no useful identity outside its parent equation, so terms
are intentionally stored as lightweight ``(node, dof, coefficient)`` tuples
rather than exposed through a separate public value-object class.  This keeps the
constraint API small while still preserving the exact user-defined term order
required by FEMaster's input syntax.

All terms, including those supplied to the constructor, pass through ``add`` so
DOF validation and numeric normalization are defined in one place.  The public
builder-style API therefore remains explicit and readable without introducing an
additional ``EquationTerm`` type solely to hold three values.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint


class Equation(Constraint):
    """Linear relation between nodal degrees of freedom."""

    def __init__(
        self,
        terms: Iterable[tuple[int | str, int, float]] = (),
    ) -> None:
        self.terms: list[tuple[int | str, int, float]] = []

        # Route constructor-provided terms through the same normalization and
        # validation path as terms appended later through the fluent API.
        for node, dof, coefficient in terms:
            self.add(node, dof, coefficient)

    def add(
        self,
        node: int | str,
        dof: int,
        coefficient: float,
    ) -> "Equation":
        """Append one ``coefficient * u(node, dof)`` term and return ``self``."""

        normalized_dof = int(dof)
        if normalized_dof < 1 or normalized_dof > 6:
            raise ValueError("Equation dof must be between 1 and 6")

        self.terms.append((node, normalized_dof, float(coefficient)))
        return self

    def export(self) -> str:
        """Export the equation using FEMaster's term-count plus triple syntax."""

        values: list[object] = []
        for node, dof, coefficient in self.terms:
            values.extend((node, dof, coefficient))

        return block([
            keyword("EQUATION"),
            csv((len(self.terms),)),
            csv(values),
        ])
