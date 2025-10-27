from __future__ import annotations

from enum import Enum
from typing import Optional, Union

from ..sets.nodeset import NodeSet

from .constraint_base import ConstraintBase, CoordinateSystemLike

class ConnectorType(str, Enum):
    NONE = "NONE"
    BEAM = "BEAM"
    HINGE = "HINGE"
    CYLINDRICAL = "CYLINDRICAL"
    TRANSLATOR = "TRANSLATOR"
    JOIN = "JOIN"
    JOINRX = "JOINRX"

    def __str__(self) -> str:  # pragma: no cover - representational
        return self.value


class Connector(ConstraintBase):
    """Connector constraint linking two reference nodes via a connector type."""

    constraint_type = "CONNECTOR"

    def __init__(
        self,
        node_set_1: NodeSet,
        node_set_2: NodeSet,
        connector_type: Union[ConnectorType, str],
        coordinate_system: Optional[CoordinateSystemLike] = None,
    ) -> None:
        super().__init__(coordinate_system)

        self.node_set_1: NodeSet = node_set_1
        self.node_set_2: NodeSet = node_set_2

        if len(self.node_set_1.nodes) != 1:
            raise ValueError("Connector node_set_1 must contain exactly one node.")
        if len(self.node_set_2.nodes) != 1:
            raise ValueError("Connector node_set_2 must contain exactly one node.")

        if not self.node_set_1.name:
            raise ValueError("Connector node_set_1 must have a name.")
        if not self.node_set_2.name:
            raise ValueError("Connector node_set_2 must have a name.")

        if isinstance(connector_type, ConnectorType):
            self.connector_type = connector_type
        else:
            value = str(connector_type).strip().upper()
            try:
                self.connector_type = ConnectorType(value)
            except ValueError as exc:  # pragma: no cover - defensive
                allowed = ", ".join(ct.value for ct in ConnectorType)
                raise ValueError(f"Unsupported connector type '{connector_type}'. Expected one of: {allowed}.") from exc

    def to_femaster(self) -> str:
        header = (
            f"*{self.constraint_type}, TYPE={self.connector_type.value}, "
            f"MASTER={self.node_set_1.name}, SLAVE={self.node_set_2.name}"
        )
        if self.coordinate_system:
            header += f", CSYS={self.coordinate_system}"
        return header

    def to_asami(self) -> str:
        raise NotImplementedError("Connector.to_asami not implemented yet.")
