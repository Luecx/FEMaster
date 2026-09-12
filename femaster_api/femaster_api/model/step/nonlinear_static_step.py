"""Nonlinear static analysis step."""

from collections.abc import Iterable

from ..common.format import block, keyword
from .constraint_method import ConstraintMethod
from .solver_control import SolverControl
from .step import Step


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
            *self.common_export_lines(),
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
