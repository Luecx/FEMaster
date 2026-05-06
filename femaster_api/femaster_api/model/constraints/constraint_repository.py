"""Repository for constraint definitions."""

from __future__ import annotations

from typing import Iterator

from .connector_constraint import ConnectorConstraint
from .coupling_constraint import CouplingConstraint
from .rbm_constraint import RBMConstraint
from .tie_constraint import TieConstraint

Constraint = RBMConstraint | CouplingConstraint | ConnectorConstraint | TieConstraint


class ConstraintRepository:
    """Repository for model-level constraints."""

    def __init__(self) -> None:
        self._items: list[Constraint] = []

    def add(self, constraint: Constraint) -> Constraint:
        if not isinstance(constraint, (RBMConstraint, CouplingConstraint, ConnectorConstraint, TieConstraint)):
            raise TypeError("constraint must be a constraint object")
        self._items.append(constraint)
        return constraint

    def all(self) -> tuple[Constraint, ...]:
        return tuple(self._items)

    def __getitem__(self, index: int | slice) -> Constraint | tuple[Constraint, ...]:
        if isinstance(index, slice):
            return tuple(self._items[index])
        return self._items[index]

    def __iter__(self) -> Iterator[Constraint]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
