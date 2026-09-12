"""Connector relation between two node regions.

A connector references two node-region names and one coordinate-system name.
Its formulation token is conceptually part of the connector definition itself,
not an independent FEMaster model object, so the known canonical tokens are
nested as ``Connector.Type`` instead of being exported through a separate
``ConnectorType`` class/module.

The constructor still accepts raw strings for forward compatibility with solver
connector formulations that may be added before the Python API is updated.  A
string matching a known canonical token is normalized to ``Connector.Type``;
unknown strings are preserved verbatim and exported unchanged.
"""

from __future__ import annotations

from enum import Enum

from ..common.format import keyword
from .constraint import Constraint


class Connector(Constraint):
    """Connector between two node regions in a named coordinate system."""

    class Type(Enum):
        """Known connector kinematic formulations supported by this API."""

        RIGID = "RIGID"
        CARTESIAN = "CARTESIAN"

    def __init__(
        self,
        type: Type | str,
        nset1: str,
        nset2: str,
        coordinate_system: str,
    ) -> None:
        # Preserve unknown native tokens while normalizing known string values to
        # the nested enum.  This keeps ``connector.type`` typed whenever possible
        # without making the API reject future FEMaster connector formulations.
        if isinstance(type, Connector.Type):
            self.type: Connector.Type | str = type
        else:
            token = str(type)
            try:
                self.type = Connector.Type(token)
            except ValueError:
                self.type = token

        self.nset1 = str(nset1)
        self.nset2 = str(nset2)
        self.coordinate_system = str(coordinate_system)

    def export(self) -> str:
        """Export this connector as one native ``*CONNECTOR`` keyword record."""

        type_token = (
            self.type.value
            if isinstance(self.type, Connector.Type)
            else self.type
        )

        return keyword(
            "CONNECTOR",
            TYPE=type_token,
            NSET1=self.nset1,
            NSET2=self.nset2,
            COORDINATESYSTEM=self.coordinate_system,
        )
