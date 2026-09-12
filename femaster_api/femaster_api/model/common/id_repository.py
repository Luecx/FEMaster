"""Ordered repository for entities identified by sparse FEMaster integer IDs.

Node and element IDs belong to the physical model and are therefore not Python
sequence positions.  ``repository[id]`` always addresses that semantic FEMaster
ID.  Positional access remains available explicitly through ``at(index)`` for
rare algorithms that need deterministic insertion order.

The generic type remains intentionally unconstrained to keep the module at one
class; stored objects are expected to expose the integer ``id`` attribute used
by node and element entities.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic, TypeVar


TId = TypeVar("TId")


class IdRepository(Generic[TId]):
    """Own sparse-ID model entities without changing their identifiers."""

    def __init__(self) -> None:
        self._items: list[TId] = []
        self._ids: dict[int, TId] = {}

    def add(self, item: TId) -> TId:
        """Insert ``item`` under its existing FEMaster ID."""

        id = item.id
        if id in self._ids:
            raise ValueError(f"duplicate id: {id}")
        self._items.append(item)
        self._ids[id] = item
        return item

    def remove(self, id: int) -> TId:
        """Remove and return one entity by semantic FEMaster ID."""

        item = self[id]
        self._items.remove(item)
        del self._ids[id]
        return item

    def __getitem__(self, id: int) -> TId:
        try:
            return self._ids[id]
        except KeyError as exc:
            raise KeyError(f"unknown id: {id}") from exc

    def __delitem__(self, id: int) -> None:
        self.remove(id)

    def at(self, index: int) -> TId:
        """Return one entity by insertion-order position."""

        try:
            return self._items[index]
        except IndexError as exc:
            raise IndexError(f"repository index out of range: {index}") from exc

    def ids(self) -> tuple[int, ...]:
        """Return semantic IDs in deterministic insertion order."""

        return tuple(item.id for item in self._items)

    def values(self) -> tuple[TId, ...]:
        """Return an immutable snapshot of owned entities."""

        return tuple(self._items)

    def __contains__(self, id: object) -> bool:
        return isinstance(id, int) and id in self._ids

    def __iter__(self) -> Iterator[TId]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
