"""Connector relation between two node regions.

A connector references two node-region names and one coordinate-system name.
The connector type accepts either the typed ``ConnectorType`` enum or a raw
native token so new solver-side connector formulations do not require a generic
Python serializer mechanism.
"""

from __future__ import annotations

from ..common.format import keyword
from .constraint import Constraint
from .constraint_connector_type import ConnectorType


class Connector(Constraint):
    """Connector between two node regions in a named coordinate system."""

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
