"""Kinematic or distributing coupling between master and slave regions.

The constraint stores semantic names only.  A slave may refer to either a node
region or a surface region; the explicit ``slave_is_surface`` flag controls the
native keyword spelling without requiring the coupling to inspect project
repositories during export.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .constraint import Constraint
from .constraint_coupling_type import CouplingType


class Coupling(Constraint):
    """Kinematic or distributing master/slave coupling."""

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
