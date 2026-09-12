"""Connector kinematic types."""

from enum import Enum


class ConnectorType(Enum):
    """Supported connector kinematic type tokens."""

    RIGID = "RIGID"
    CARTESIAN = "CARTESIAN"
