"""Linear static structural analysis step.

``StaticStep`` extends the common loadcase definition with the optional inertia
relief and load-rebalancing controls implemented by FEMaster.  Load and support
collectors, solver selection and constraint treatment are inherited from
``Step``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, keyword
from .step import Step
from .util.constraint_method import ConstraintMethod
from .util.solver_control import SolverControl


class StaticStep(Step):
    """Linear static structural analysis."""

    type_name = "LINEARSTATIC"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
        inertia_relief: bool = False,
        rebalance_loads: bool = False,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
            constraint_method=constraint_method,
        )
        self.inertia_relief = bool(inertia_relief)
        self.rebalance_loads = bool(rebalance_loads)

    def export(self) -> str:
        lines = self._common_lines()
        if self.inertia_relief:
            lines.append(keyword("INERTIARELIEF"))
        if self.rebalance_loads:
            lines.append(keyword("REBALANCELOADS"))
        return block(lines)
