"""Reusable FEMaster parts and their protected repository.

``Part`` owns all local topology and assignments.  ``PartRepository`` additionally
represents the implicit root/default part permanently at index zero, preventing a
second copy of default-part state on ``Project`` and preserving one consistent
ownership path.
"""

from .part import Part
from .part_repository import PartRepository

__all__ = ["Part", "PartRepository"]
