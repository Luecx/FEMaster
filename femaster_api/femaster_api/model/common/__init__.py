"""Small infrastructure shared by independent FEMaster model concepts.

Only generic ownership, identifier and text-format helpers live here.  Domain
relationships are deliberately expressed with the concrete model classes at the
use site rather than through generic ``EntityReference`` / ``NodeReference`` /
``ElementReference`` aliases.  This keeps every public signature explicit about
which FEMaster objects are actually legal.
"""

from .id_repository import IdRepository
from .named_object import NamedObject
from .named_repository import NamedRepository

__all__ = [
    "IdRepository",
    "NamedObject",
    "NamedRepository",
]
