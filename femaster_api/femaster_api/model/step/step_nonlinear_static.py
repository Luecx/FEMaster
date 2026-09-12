"""Nonlinear static structural analysis with concrete collector references.

The class exposes FEMaster's nonlinear increment controls directly.  Loads and
supports are object relationships inherited from ``Step``; only the nonlinear
control mode itself remains a native token because it is intrinsic procedure
configuration rather than a reference to another model object.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, keyword
from ..load.load_collector import LoadCollector
from ..support.support_collector import SupportCollector
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
        loads: Iterable[LoadCollector] = (),
        supports: Iterable[SupportCollector] = (),
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
