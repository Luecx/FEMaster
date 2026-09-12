"""Ordered heterogeneous repository of assembly-level constraints.

Constraints generally have no shared semantic name, so a list-like repository is
more appropriate than ``NamedRepository``.  The repository preserves definition
order and delegates export to each concrete constraint without inspecting its
type.
"""

from __future__ import annotations

from collections.abc import Iterator

from .constraint import Constraint


class ConstraintRepository:
    """Own assembly constraints in deterministic insertion order."""

    def __init__(self) -> None:
        self._items: list[Constraint] = []

    def add(self, constraint: Constraint) -> Constraint:
        """Append a concrete constraint and return it."""

        if not isinstance(constraint, Constraint):
            raise TypeError("constraint must derive from Constraint")
        self._items.append(constraint)
        return constraint

    def export(self) -> str:
        return "\n\n".join(item.export() for item in self._items)

    def __getitem__(self, index: int) -> Constraint:
        return self._items[index]

    def __iter__(self) -> Iterator[Constraint]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
