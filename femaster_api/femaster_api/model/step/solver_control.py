"""Shared linear solver control."""

from ..common.format import keyword
from .solver_device import SolverDevice
from .solver_method import SolverMethod


class SolverControl:
    """Linear solver device/method selection."""

    def __init__(
        self,
        device: SolverDevice = SolverDevice.CPU,
        method: SolverMethod = SolverMethod.DIRECT,
    ) -> None:
        self.device = device
        self.method = method

    def export(self) -> str:
        return keyword("SOLVER", DEVICE=self.device.value, METHOD=self.method.value)
