"""Canonical connector kinematic type tokens.

FEMaster may grow additional connector types; the enum contains the currently
modeled common forms while ``Connector`` also accepts a raw string for native
extensions not yet represented here.
"""

from __future__ import annotations

from enum import Enum


class ConnectorType(Enum):
    """Known connector kinematic formulations."""

    RIGID = "RIGID"
    CARTESIAN = "CARTESIAN"
