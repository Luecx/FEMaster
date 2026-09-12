"""Base object for model entities with immutable semantic names.

Named model objects are indexed by ``NamedRepository``.  Their name therefore
forms part of repository identity and must remain stable after insertion.
Making only the name immutable keeps the surrounding FEM definition mutable
without allowing repository lookup tables to become stale.
"""

from __future__ import annotations


class NamedObject:
    """Base class for FEMaster objects identified by a stable semantic name."""

    __slots__ = ("_name",)

    def __init__(self, name: str) -> None:
        normalized = str(name).strip()
        if not normalized:
            raise ValueError("name must not be empty")
        self._name = normalized

    @property
    def name(self) -> str:
        """Return the immutable semantic name used for repository lookup."""

        return self._name

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r})"
