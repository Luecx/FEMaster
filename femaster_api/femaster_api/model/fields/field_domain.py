"""Field domain enum."""

from __future__ import annotations

from enum import Enum


class FieldDomain(Enum):
    UNKNOWN = "UNKNOWN"
    NODE = "NODE"
    ELEMENT = "ELEMENT"
    ELEMENT_NODAL = "ELEMENT_NODAL"
    ELEMENT_IP = "ELEMENT_IP"
