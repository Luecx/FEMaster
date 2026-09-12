"""Linear buckling analysis step."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .solver_control import SolverControl
from .step import Step


class BucklingStep(Step):
    """Linearized eigenvalue buckling analysis."""

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
        super().__init__(name, loads=loads, supports=supports, solver=solver)
        self.number_of_modes = int(number_of_modes)
        self.sigma = None if sigma is None else float(sigma)

    def export(self) -> str:
        lines = [
            *self.common_export_lines(),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ]
        if self.sigma is not None:
            lines.extend([keyword("SIGMA"), csv((self.sigma,))])
        return block(lines)
