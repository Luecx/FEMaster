"""Linear static topology loadcase definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.fields import Field
from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class TopologyStaticLoadcase:
    """Linear static topology optimization loadcase."""

    name: str
    density: Field
    orientation: Field | None = None
    exponent: float | None = None
    loads: tuple[LoadCollector, ...] = ()
    supports: tuple[SupportCollector, ...] = ()
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE

    def __post_init__(self) -> None:
        object.__setattr__(self, "loads", tuple(self.loads))
        object.__setattr__(self, "supports", tuple(self.supports))
