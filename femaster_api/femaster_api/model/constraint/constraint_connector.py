"""Connector relation between concrete node regions and a coordinate system.

A connector stores two ``NodeRegion`` objects and one ``CoordinateSystem``
object.  Their semantic names are native-format details and are never accepted as
constructor substitutes.  The kinematic formulation remains the nested
``Connector.Type`` enum because it belongs to the connector definition rather
than to the global model graph.
"""

from __future__ import annotations

from enum import Enum

from ..common.format import keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..region.region_node import NodeRegion
from .constraint import Constraint


class Connector(Constraint):
    """Connector between two node regions in one coordinate system."""

    class Type(Enum):
        """Known connector kinematic formulations supported by this API."""

        RIGID = "RIGID"
        CARTESIAN = "CARTESIAN"

    def __init__(
        self,
        type: Type,
        nset1: NodeRegion,
        nset2: NodeRegion,
        coordinate_system: CoordinateSystem,
    ) -> None:
        if not isinstance(type, Connector.Type):
            raise TypeError("type must be Connector.Type")
        if not isinstance(nset1, NodeRegion) or not isinstance(nset2, NodeRegion):
            raise TypeError("nset1 and nset2 must be NodeRegion objects")
        if not isinstance(coordinate_system, CoordinateSystem):
            raise TypeError("coordinate_system must be a CoordinateSystem object")

        self.type = type
        self.nset1 = nset1
        self.nset2 = nset2
        self.coordinate_system = coordinate_system

    def export(self) -> str:
        """Export this connector as one native ``*CONNECTOR`` keyword record."""

        return keyword(
            "CONNECTOR",
            TYPE=self.type.value,
            NSET1=self.nset1.name,
            NSET2=self.nset2.name,
            COORDINATESYSTEM=self.coordinate_system.name,
        )
