"""Solver control data object."""

from __future__ import annotations

from dataclasses import dataclass

from .solver_device import SolverDevice
from .solver_method import SolverMethod


@dataclass(frozen=True, slots=True)
class SolverControl:
    device: SolverDevice = SolverDevice.CPU
    method: SolverMethod = SolverMethod.DIRECT
