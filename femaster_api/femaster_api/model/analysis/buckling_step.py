"""Linear buckling analysis step."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class BucklingStep:
    name: str
    loads: tuple[LoadCollector, ...]
    supports: tuple[SupportCollector, ...] = ()
    number_of_modes: int = 10
    sigma: float | None = None
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
    request_stiffness: str | None = None
    request_stgeom: str | None = None
