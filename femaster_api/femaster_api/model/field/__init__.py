"""Central field representation shared by model input and result data.

``FieldDomain`` defines physical storage location, while ``FieldType`` defines
semantic meaning independently of file format.  Editable project fields and
post-processing result fields therefore use the same core representation and do
not maintain competing RES/FRD-specific enums.
"""

from .field import Field
from .field_domain import FieldDomain
from .field_repository import FieldRepository
from .field_type import FieldType
from .typing import FieldIndex, FieldKey

__all__ = [
    "Field",
    "FieldDomain",
    "FieldIndex",
    "FieldKey",
    "FieldRepository",
    "FieldType",
]
