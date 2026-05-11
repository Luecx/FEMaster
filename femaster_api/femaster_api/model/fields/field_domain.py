"""Field domain enum."""

from __future__ import annotations

from enum import Enum


class FieldDomain(Enum):
    UNKNOWN = "UNKNOWN"
    NODE = "NODE"
    ELEMENT = "ELEMENT"
    ELEMENT_NODAL = "ELEMENTNODAL"
    IP = "IP"
