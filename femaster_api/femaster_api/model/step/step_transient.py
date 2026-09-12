"""Linear structural transient step with concrete load/support collectors.

The step owns physical time control and optionally Newmark parameters, Rayleigh
damping and result-write cadence.  Load and support relationships are actual
collector objects inherited from ``Step``; their names appear only in exported
native input syntax.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..load.load_collector import LoadCollector
from ..support.support_collector import SupportCollector
from .step import Step
from .util.newmark_control import NewmarkControl
from .util.rayleigh_damping import RayleighDamping
from .util.solver_control import SolverControl
from .util.time_control import TimeControl


class TransientStep(Step):
    """Linear transient structural analysis."""

    type_name = "LINEARTRANSIENT"

    def __init__(
        self,
        name: str,
        time: TimeControl,
        *,
        loads: Iterable[LoadCollector] = (),
        supports: Iterable[SupportCollector] = (),
        solver: SolverControl | None = None,
        newmark: NewmarkControl | None = None,
        damping: RayleighDamping | None = None,
        write_every: int | None = None,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
        )
        self.time = time
        self.newmark = newmark
        self.damping = damping
        self.write_every = write_every

    def export(self) -> str:
        lines = [
            *self._common_lines(),
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
