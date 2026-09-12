"""Nonlinear static structural analysis step.

The class exposes FEMaster's nonlinear increment controls directly instead of
hiding them in an opaque options mapping.  ``None`` means that the solver default
is retained.  The control mode is stored as the native upper-case token so load
control and arc-length-style extensions can remain visible in the exported deck.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, keyword
from .step import Step
from .util.constraint_method import ConstraintMethod
from .util.solver_control import SolverControl


class NonlinearStaticStep(Step):
    """Geometrically/materially nonlinear static analysis."""

    type_name = "NONLINEARSTATIC"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
        control: str = "LOAD",
        increments: int | None = None,
        max_increments: int | None = None,
        initial_increment: float | None = None,
        minimum_increment: float | None = None,
        maximum_increment: float | None = None,
        max_iterations: int | None = None,
        tolerance: float | None = None,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
            constraint_method=constraint_method,
        )
        self.control = str(control).upper()
        self.increments = increments
        self.max_increments = max_increments
        self.initial_increment = initial_increment
        self.minimum_increment = minimum_increment
        self.maximum_increment = maximum_increment
        self.max_iterations = max_iterations
        self.tolerance = tolerance

    def export(self) -> str:
        return block([
            *self._common_lines(),
            keyword(
                "NONLINEAR",
                CONTROL=self.control,
                INCREMENTS=self.increments,
                MAX_INCREMENTS=self.max_increments,
                INITIAL_INCREMENT=self.initial_increment,
                MINIMUM_INCREMENT=self.minimum_increment,
                MAXIMUM_INCREMENT=self.maximum_increment,
                MAXITER=self.max_iterations,
                TOL=self.tolerance,
            ),
        ])
