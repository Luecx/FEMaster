from __future__ import annotations

from typing import Iterable, Optional, Tuple

from sets.nodeset import NodeSet

from .constraint_base import ConstraintBase

class Coupling(ConstraintBase):
    """Base representation for couplings relating a master to a slave nodeset."""

    coupling_type: str = "COUPLING"

    def __init__(self, master: NodeSet, slave: NodeSet) -> None:
        super().__init__(None)
        if not isinstance(master, NodeSet):
            raise TypeError("Coupling master must be provided as a NodeSet.")
        if not isinstance(slave, NodeSet):
            raise TypeError("Coupling slave must be provided as a NodeSet.")

        if not master.name:
            raise ValueError("Coupling master nodeset must have a name.")
        if not slave.name:
            raise ValueError("Coupling slave nodeset must have a name.")

        if len(master.nodes) != 1:
            raise ValueError("Coupling master nodeset must contain exactly one node.")
        if not slave.nodes:
            raise ValueError("Coupling slave nodeset must contain at least one node.")

        self.master = master
        self.slave = slave

    @staticmethod
    def _token(nodeset: NodeSet) -> str:
        if not nodeset.name:
            raise ValueError("Coupling nodeset must have a name for export.")
        return nodeset.name


class KinematicCoupling(Coupling):
    """Kinematic coupling with explicit DOF activation flags."""

    coupling_type = "KINEMATIC"

    def __init__(self, master: NodeSet, slave: NodeSet, dofs: Iterable[int | bool]) -> None:
        super().__init__(master, slave)

        vals = tuple(1 if bool(v) else 0 for v in dofs)
        if len(vals) != 6:
            raise ValueError("Kinematic coupling requires 6 DOF flags (ux, uy, uz, rx, ry, rz).")
        self.dofs: Tuple[int, int, int, int, int, int] = vals

    def to_femaster(self) -> str:
        master = self._token(self.master)
        slave = self._token(self.slave)
        header = f"*COUPLING, TYPE={self.coupling_type}, MASTER={master}, SLAVE={slave}"
        payload = ", ".join(str(v) for v in self.dofs)
        return "\n".join((header, payload))

    def to_asami(self) -> str:
        raise NotImplementedError("KinematicCoupling.to_asami not implemented yet.")


class DistributingCoupling(Coupling):
    """Distributing coupling optionally carrying weighting information."""

    coupling_type = "DISTRIBUTING"

    def __init__(
        self,
        master: NodeSet,
        slave: NodeSet,
        weighting: Optional[str] = None,
    ) -> None:
        super().__init__(master, slave)
        self.weighting: Optional[str] = (weighting.strip().upper() if weighting else None)

    def to_femaster(self) -> str:
        master = self._token(self.master)
        slave = self._token(self.slave)
        header = f"*COUPLING, TYPE={self.coupling_type}, MASTER={master}, SLAVE={slave}"
        if self.weighting:
            header += f", WEIGHTING={self.weighting}"
        return header

    def to_asami(self) -> str:
        raise NotImplementedError("DistributingCoupling.to_asami not implemented yet.")
