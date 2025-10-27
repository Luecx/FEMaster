from __future__ import annotations

from typing import Optional

from ..sets.nodeset import NodeSet

from .constraint_base import ConstraintBase

class Tie(ConstraintBase):
    """Tie constraint linking a master reference node to a slave nodeset."""

    constraint_type = "TIE"

    def __init__(
        self,
        master: NodeSet,
        slave: NodeSet,
        *,
        adjust: Optional[bool] = None,
        distance: Optional[float] = None,
    ) -> None:
        super().__init__(None)

        if not isinstance(master, NodeSet):
            raise TypeError("Tie master must be provided as a NodeSet.")
        if not isinstance(slave, NodeSet):
            raise TypeError("Tie slave must be provided as a NodeSet.")

        if not master.name:
            raise ValueError("Tie master nodeset must have a name.")
        if not slave.name:
            raise ValueError("Tie slave nodeset must have a name.")

        if len(master.nodes) != 1:
            raise ValueError("Tie master nodeset must contain exactly one node.")
        if not slave.nodes:
            raise ValueError("Tie slave nodeset must contain at least one node.")

        self.master = master
        self.slave = slave

        self.adjust: Optional[bool] = (None if adjust is None else bool(adjust))
        self.distance: Optional[float] = (None if distance is None else float(distance))

    @staticmethod
    def _token(nodeset: NodeSet) -> str:
        if not nodeset.name:
            raise ValueError("Tie nodeset must have a name for export.")
        return nodeset.name

    def to_femaster(self) -> str:
        master = self._token(self.master)
        slave = self._token(self.slave)

        header = f"*{self.constraint_type}, MASTER={master}, SLAVE={slave}"
        if self.distance is not None:
            header += f", DISTANCE={self.distance}"
        if self.adjust is not None:
            header += f", ADJUST={'YES' if self.adjust else 'NO'}"
        return header

    def to_asami(self) -> str:
        raise NotImplementedError("Tie.to_asami not implemented yet.")
