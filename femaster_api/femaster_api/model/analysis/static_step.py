"""Linear static analysis step."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class StaticStep:
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
