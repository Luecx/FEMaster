"""Shared solver-selection control for structural analysis steps.

``SolverControl`` combines the orthogonal device and method selections and owns
the corresponding native ``*SOLVER`` keyword.  Keeping the control outside the
concrete step classes avoids repeating identical export logic while still
representing a genuine shared FEMaster concept.
"""

from __future__ import annotations

from ...common.format import keyword
from .solver_device import SolverDevice
from .solver_method import SolverMethod


class SolverControl:
    """Device/method pair used to solve one analysis system."""

    def __init__(
        self,
        device: SolverDevice = SolverDevice.CPU,
        method: SolverMethod = SolverMethod.DIRECT,
    ) -> None:
        self.device = device
        self.method = method

    def export(self) -> str:
        """Export the native solver-selection keyword."""

        return keyword(
            "SOLVER",
            DEVICE=self.device.value,
            METHOD=self.method.value,
        )
