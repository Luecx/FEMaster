"""Connector type enum."""

from __future__ import annotations

from enum import Enum


class ConnectorType(Enum):
    BEAM = "BEAM"
    HINGE = "HINGE"
    CYLINDRICAL = "CYLINDRICAL"
    TRANSLATOR = "TRANSLATOR"
    JOIN = "JOIN"
    JOINRX = "JOINRX"
