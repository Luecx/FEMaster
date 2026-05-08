"""Linear static loadcase definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class LinearStaticLoadcase:
    """Linear static analysis loadcase."""

    name: str
    loads: tuple[LoadCollector, ...] = ()
    supports: tuple[SupportCollector, ...] = ()
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
    inertia_relief: bool = False
    inertia_relief_consider_point_masses: bool = True
    rebalance_loads: bool = False
    request_stiffness: str | None = None
    constraint_summary: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "loads", tuple(self.loads))
        object.__setattr__(self, "supports", tuple(self.supports))
