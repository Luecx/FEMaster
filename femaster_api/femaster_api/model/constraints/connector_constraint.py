"""Connector constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import NodeSet

from .connector_type import ConnectorType


@dataclass(frozen=True, slots=True)
class ConnectorConstraint:
    """Connector between two node sets."""

    type: ConnectorType
    nset1: NodeSet
    nset2: NodeSet
    coordinate_system: Orientation
