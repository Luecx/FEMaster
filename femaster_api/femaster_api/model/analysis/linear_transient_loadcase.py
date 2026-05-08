"""Linear transient loadcase definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.fields import Field
from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .newmark_control import NewmarkControl
from .rayleigh_damping import RayleighDamping
from .solver_control import SolverControl
from .time_control import TimeControl


@dataclass(frozen=True, slots=True)
class LinearTransientLoadcase:
    """Linear transient analysis loadcase."""

    name: str
    loads: tuple[LoadCollector, ...]
    time: TimeControl
    supports: tuple[SupportCollector, ...] = ()
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
    newmark: NewmarkControl | None = None
    damping: RayleighDamping | None = None
    write_every_type: str | None = None
    write_every: float | None = None
    initial_velocity: Field | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "loads", tuple(self.loads))
        object.__setattr__(self, "supports", tuple(self.supports))
