"""Kinematic/distributing coupling."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint
from .coupling_type import CouplingType


class Coupling(Constraint):
    """Coupling between master and slave node/surface regions."""

    def __init__(
        self,
        master: str,
        slave: str,
        *,
        type: CouplingType = CouplingType.KINEMATIC,
        dofs: Iterable[bool | int] = (1, 1, 1, 1, 1, 1),
        slave_is_surface: bool = False,
    ) -> None:
        self.master = str(master)
        self.slave = str(slave)
        self.type = type
        self.dofs = tuple(int(bool(value)) for value in dofs)
        self.slave_is_surface = bool(slave_is_surface)
        if len(self.dofs) != 6:
            raise ValueError("Coupling requires exactly 6 DOF flags")

    def export(self) -> str:
        return block([
            keyword(
                "COUPLING",
                MASTER=self.master,
                TYPE=self.type.value,
                SFSET=self.slave if self.slave_is_surface else None,
                SLAVE=None if self.slave_is_surface else self.slave,
            ),
            csv(self.dofs),
        ])
