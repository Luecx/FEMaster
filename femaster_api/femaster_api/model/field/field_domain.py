"""Field storage domains."""

from enum import Enum


class FieldDomain(Enum):
    """Physical storage domains implemented by FEMaster ModelData."""

    UNKNOWN = "UNKNOWN"
    NODE = "NODE"
    ELEMENT = "ELEMENT"
    ELEMENT_NODAL = "ELEMENT_NODAL"
    ELEMENT_IP = "ELEMENT_IP"
    ELEMENT_MP = "ELEMENT_MP"
