"""Entity set domain enum."""

from __future__ import annotations

from enum import Enum


class EntityType(Enum):
    NODE = "node"
    ELEMENT = "element"
    SURFACE = "surface"
