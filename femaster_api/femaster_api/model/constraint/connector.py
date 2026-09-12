"""Connector relation."""

from ..common.format import keyword
from .connector_type import ConnectorType
from .constraint import Constraint


class Connector(Constraint):
    """Connector relation between two node regions."""

    def __init__(
        self,
        type: ConnectorType | str,
        nset1: str,
        nset2: str,
        coordinate_system: str,
    ) -> None:
        self.type = type.value if isinstance(type, ConnectorType) else str(type)
        self.nset1 = str(nset1)
        self.nset2 = str(nset2)
        self.coordinate_system = str(coordinate_system)

    def export(self) -> str:
        return keyword(
            "CONNECTOR",
            TYPE=self.type,
            NSET1=self.nset1,
            NSET2=self.nset2,
            COORDINATESYSTEM=self.coordinate_system,
        )
