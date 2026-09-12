"""Small infrastructure shared by independent FEMaster model concepts.

Only generic ownership, identifier and text-format helpers live here.  No
finite-element concept such as nodes, elements, regions, steps or results should
be introduced in ``common`` merely to avoid choosing its proper domain folder.
"""

from .id_repository import IdRepository
from .named_object import NamedObject
from .named_repository import NamedRepository
from .typing import ElementReference, EntityReference, NodeReference

__all__ = [
    "ElementReference",
    "EntityReference",
    "IdRepository",
    "NamedObject",
    "NamedRepository",
    "NodeReference",
]
