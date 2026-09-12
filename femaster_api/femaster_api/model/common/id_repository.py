"""Ordered repository preserving sparse FEMaster integer IDs."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic, TypeVar


T = TypeVar("T")


class IdRepository(Generic[T]):
    """Repository whose integer subscription addresses the semantic FEMaster ID."""

    def __init__(self) -> None:
        self._items: list[T] = []
        self._ids: dict[int, T] = {}

    def add(self, item: T) -> T:
        if item.id in self._ids:
            raise ValueError(f"duplicate id: {item.id}")
        self._items.append(item)
        self._ids[item.id] = item
        return item

    def remove(self, id: int) -> T:
        item = self[id]
        self._items.remove(item)
        del self._ids[id]
        return item

    def __getitem__(self, id: int) -> T:
        try:
            return self._ids[id]
        except KeyError as exc:
            raise KeyError(f"unknown id: {id}") from exc

    def __delitem__(self, id: int) -> None:
        self.remove(id)

    def at(self, index: int) -> T:
        try:
            return self._items[index]
        except IndexError as exc:
            raise IndexError(f"repository index out of range: {index}") from exc

    def ids(self) -> tuple[int, ...]:
        return tuple(item.id for item in self._items)

    def values(self) -> tuple[T, ...]:
        return tuple(self._items)

    def __contains__(self, id: object) -> bool:
        return isinstance(id, int) and id in self._ids

    def __iter__(self) -> Iterator[T]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
