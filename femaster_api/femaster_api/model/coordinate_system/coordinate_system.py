"""Base class for named global coordinate-system definitions.

Coordinate systems are project-level reusable objects used directly by loads,
supports, sections and constraints.  Consumers store the ``CoordinateSystem``
object itself; the immutable semantic name is only a native serialization token.
Concrete subclasses own the geometry needed by one orientation type and their
exact export syntax.
"""

from __future__ import annotations

from ..common.named_object import NamedObject


class CoordinateSystem(NamedObject):
    """Base class for one named FEMaster coordinate system."""

    def export(self) -> str:
        raise NotImplementedError
