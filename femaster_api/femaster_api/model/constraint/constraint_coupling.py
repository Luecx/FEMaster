"""Kinematic/distributing coupling between concrete master and slave objects.

The master is a ``Node`` or ``NodeRegion``.  The slave is a ``Node`` /
``NodeRegion`` or a ``Surface`` / ``SurfaceRegion``.  The native distinction
between ``SLAVE=`` and ``SFSET=`` is inferred from the actual Python type rather
than stored in a parallel boolean flag or encoded in a string reference.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..node.node import Node
from ..region.region_node import NodeRegion
from ..region.region_surface import SurfaceRegion
from ..surface.surface import Surface
from .constraint import Constraint
from .constraint_coupling_type import CouplingType


class Coupling(Constraint):
    """Kinematic or distributing master/slave coupling."""

    def __init__(
        self,
        master: Node | NodeRegion,
        slave: Node | NodeRegion | Surface | SurfaceRegion,
        *,
        type: CouplingType = CouplingType.KINEMATIC,
        dofs: Iterable[bool | int] = (1, 1, 1, 1, 1, 1),
    ) -> None:
        if not isinstance(master, (Node, NodeRegion)):
            raise TypeError("master must be Node or NodeRegion")
        if not isinstance(slave, (Node, NodeRegion, Surface, SurfaceRegion)):
            raise TypeError(
                "slave must be Node, NodeRegion, Surface or SurfaceRegion"
            )
        if not isinstance(type, CouplingType):
            raise TypeError("type must be CouplingType")

        self.master = master
        self.slave = slave
        self.type = type
        self.dofs = tuple(int(bool(value)) for value in dofs)

        if len(self.dofs) != 6:
            raise ValueError("Coupling requires exactly 6 DOF flags")

    def export(self) -> str:
        master = self.master.id if isinstance(self.master, Node) else self.master.name
        slave_is_surface = isinstance(self.slave, (Surface, SurfaceRegion))
        slave = (
            self.slave.id
            if isinstance(self.slave, Node)
            else self.slave.name
        )

        return block([
            keyword(
                "COUPLING",
                MASTER=master,
                TYPE=self.type.value,
                SFSET=slave if slave_is_surface else None,
                SLAVE=None if slave_is_surface else slave,
            ),
            csv(self.dofs),
        ])
