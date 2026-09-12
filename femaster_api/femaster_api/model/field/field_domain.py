"""Physical storage domains shared by model fields and result fields.

The enum mirrors the domains implemented by FEMaster ``ModelData``.  It is the
single domain definition used by input fields, native RES parsing and FRD result
conversion, preventing each file format from inventing a competing set of
location enums.
"""

from __future__ import annotations

from enum import Enum


class FieldDomain(Enum):
    """Location at which each row of a field is defined."""

    UNKNOWN = "UNKNOWN"
    NODE = "NODE"
    ELEMENT = "ELEMENT"
    ELEMENT_NODAL = "ELEMENT_NODAL"
    ELEMENT_IP = "ELEMENT_IP"
    ELEMENT_MP = "ELEMENT_MP"
