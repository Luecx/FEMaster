"""Nonlinear static analysis step."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.loads import LoadCollector
from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


@dataclass(frozen=True, slots=True)
class NonlinearStaticStep:
    name: str
    loads: tuple[LoadCollector, ...] = ()
    supports: tuple[SupportCollector, ...] = ()
    solver: SolverControl | None = None
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
    control: str = "LOAD"
    increments: int | None = None
    max_increments: int | None = None
    initial_increment: float | None = None
    minimum_increment: float | None = None
    maximum_increment: float | None = None
    arc_length_psi: float | None = None
    adaptive: bool | None = None
    growth_factor: float | None = None
    cutback_factor: float | None = None
    fast_iterations: int | None = None
    slow_iterations: int | None = None
    maximum_cutbacks: int | None = None
    max_iterations: int | None = None
    tolerance: float | None = None
    regularize_zero_rows: bool | None = None
    regularization_alpha: float | None = None
    request_stiffness: str | None = None
    constraint_summary: bool = False
