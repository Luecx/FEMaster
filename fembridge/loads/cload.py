from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from ..nodes.node import Node
from ..sets.nodeset import NodeSet

from .load_base import Load, CoordinateSystem, AmplitudeFunction, TIME_DEPENDENT_ERROR

Target = Union[Node, NodeSet]


class CLoad(Load):
    """Concentrated load on either a single Node or a NodeSet."""

    load_type = "CLOAD"

    def __init__(
        self,
        target: Target,
        comps: Iterable[float],
        coordinate_system: Optional[CoordinateSystem] = None,
        amplitude: Optional[AmplitudeFunction] = None,
    ) -> None:
        super().__init__(coordinate_system=coordinate_system, amplitude=amplitude)
        self.nodeset: NodeSet = NodeSet.internal(target)

        values = tuple(float(x) for x in comps)
        if len(values) not in (3, 6):
            raise ValueError("CLoad expects 3 (Fx,Fy,Fz) or 6 (Fx,Fy,Fz,Mx,My,Mz) components.")
        if len(values) == 3:
            values = (*values, 0.0, 0.0, 0.0)
        self.components: Tuple[float, float, float, float, float, float] = values

    def to_femaster(self, collector_name: str) -> str:
        if self.is_time_dependent:
            raise NotImplementedError(TIME_DEPENDENT_ERROR)
        fx, fy, fz, mx, my, mz = self.components
        token   = self.nodeset.name or str(int(self.nodeset.nodes[0].node_id))
        payload = f"{token}, {fx}, {fy}, {fz}, {mx}, {my}, {mz}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("CLoad.to_asami not implemented yet.")
