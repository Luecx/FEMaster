from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from ..nodes.node import Node
from ..sets.nodeset import NodeSet

from .support_base import SupportBase, CoordinateSystemLike


Target = Union[Node, NodeSet]


class Support(SupportBase):
    """Support entry applied to a Node or NodeSet with optional coordinate system."""

    support_type = "SUPPORT"

    def __init__(
        self,
        target: Target,
        dofs: Iterable[Optional[int]],
        coordinate_system: Optional[CoordinateSystemLike] = None,
    ) -> None:
        super().__init__(coordinate_system)

        vals = tuple(v if v is None else int(v) for v in dofs)
        if len(vals) != 6:
            raise ValueError("dofs must have length 6 (UX, UY, UZ, RX, RY, RZ).")
        for v in vals:
            if v is not None and v not in (0, 1):
                raise ValueError("Each dof must be 0, 1, or None.")

        self.nodeset: NodeSet = NodeSet.internal(target)
        self.dofs: Tuple[Optional[int], Optional[int], Optional[int], Optional[int], Optional[int], Optional[int]] = vals

    def to_femaster(self, collector_name: str) -> str:
        fields = [("" if v is None else str(int(v))) for v in self.dofs]
        token = self.nodeset.name or str(int(self.nodeset.nodes[0].node_id))
        payload = f"{token}, " + ", ".join(fields)
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("Support.to_asami not implemented yet.")
