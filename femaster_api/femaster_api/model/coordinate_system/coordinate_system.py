"""Base class for named global coordinate-system definitions.

Coordinate systems are project-level reusable definitions referenced by loads,
supports, sections and other model concepts through semantic names.  Concrete
subclasses own the geometry needed by one native orientation type and therefore
also own their exact export syntax.
"""

from __future__ import annotations

from ..common.named_object import NamedObject


class CoordinateSystem(NamedObject):
    """Base class for one named FEMaster coordinate system."""

    def export(self) -> str:
        raise NotImplementedError
