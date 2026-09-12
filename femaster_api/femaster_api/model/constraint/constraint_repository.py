"""Ordered heterogeneous constraint repository."""

from collections.abc import Iterator

from .constraint import Constraint


class ConstraintRepository:
    """Ordered collection of heterogeneous assembly constraints."""

    def __init__(self) -> None:
        self._items: list[Constraint] = []

    def add(self, constraint: Constraint) -> Constraint:
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
