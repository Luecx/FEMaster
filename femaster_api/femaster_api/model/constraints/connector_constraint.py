"""Connector constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import EntitySet

from .connector_type import ConnectorType


@dataclass(frozen=True, slots=True)
class ConnectorConstraint:
    type: ConnectorType
    nset1: EntitySet
    nset2: EntitySet
    coordinate_system: Orientation
