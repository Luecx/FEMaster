"""Linear buckling loadcase definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class LinearBucklingLoadcase:
    """Linear buckling analysis loadcase."""

    name: str
    loads: tuple[LoadCollector, ...]
    supports: tuple[SupportCollector, ...] = ()
    number_of_modes: int = 10
    sigma: float | None = None
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
    request_stiffness: str | None = None
    request_stgeom: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "loads", tuple(self.loads))
        object.__setattr__(self, "supports", tuple(self.supports))
