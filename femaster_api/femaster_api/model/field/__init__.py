"""Central field representation shared by model input and result data.

``FieldDomain`` defines physical storage location, while ``FieldType`` defines
semantic meaning independently of file format.  Editable model fields use
object-valued entity addresses; result readers preserve serialized solver
addresses internally because results do not own the editable model graph.
"""

from .field import Field
from .field_domain import FieldDomain
from .field_repository import FieldRepository
from .field_type import FieldType

__all__ = [
    "Field",
    "FieldDomain",
    "FieldRepository",
    "FieldType",
]
