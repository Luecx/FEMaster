"""Coordinate-system type enum."""

from __future__ import annotations

from enum import Enum


class OrientationType(Enum):
    RECTANGULAR = "RECTANGULAR"
    CYLINDRICAL = "CYLINDRICAL"
