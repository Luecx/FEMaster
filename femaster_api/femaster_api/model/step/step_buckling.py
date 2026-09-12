"""Linearized eigenvalue buckling analysis step.

The procedure combines a static preload definition with a requested number of
buckling eigenpairs.  An optional spectral shift is emitted through ``*SIGMA``.
All collector and solver references remain semantic names inherited from
``Step``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .step import Step
from .util.solver_control import SolverControl


class BucklingStep(Step):
    """Linearized buckling eigenvalue analysis."""

    type_name = "LINEARBUCKLING"

    def __init__(
        self,
        name: str,
        number_of_modes: int,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        sigma: float | None = None,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
        )
        self.number_of_modes = int(number_of_modes)
        self.sigma = None if sigma is None else float(sigma)

    def export(self) -> str:
        lines = [
            *self._common_lines(),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ]
        if self.sigma is not None:
            lines.extend([keyword("SIGMA"), csv((self.sigma,))])
        return block(lines)
