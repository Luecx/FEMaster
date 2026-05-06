"""Field domain enum."""

from __future__ import annotations

from enum import Enum


class FieldDomain(Enum):
    NODE = "NODE"
    ELEMENT = "ELEMENT"
    IP = "IP"
