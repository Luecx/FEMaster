"""Linear transient analysis step."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .newmark_control import NewmarkControl
from .rayleigh_damping import RayleighDamping
from .solver_control import SolverControl
from .step import Step
from .time_control import TimeControl


class TransientStep(Step):
    """Linear structural transient analysis with Newmark integration."""

    type_name = "LINEARTRANSIENT"

    def __init__(
        self,
        name: str,
        time: TimeControl,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        newmark: NewmarkControl | None = None,
        damping: RayleighDamping | None = None,
        write_every: int | None = None,
    ) -> None:
        super().__init__(name, loads=loads, supports=supports, solver=solver)
        self.time = time
        self.newmark = newmark
        self.damping = damping
        self.write_every = write_every

    def export(self) -> str:
        lines = [
            *self.common_export_lines(),
            keyword("TIME"),
            csv((self.time.start, self.time.end, self.time.step)),
        ]
        if self.newmark is not None:
            lines.extend([
                keyword("NEWMARK"),
                csv((self.newmark.beta, self.newmark.gamma)),
            ])
        if self.damping is not None:
            lines.extend([
                keyword("DAMPING", TYPE="RAYLEIGH"),
                csv((self.damping.alpha, self.damping.beta)),
            ])
        if self.write_every is not None:
            lines.extend([
                keyword("WRITEEVERY", TYPE="STEPS"),
                csv((self.write_every,)),
            ])
        return block(lines)
