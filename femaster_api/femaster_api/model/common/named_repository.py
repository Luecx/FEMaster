"""Ordered name- and index-addressable repository."""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic, TypeVar, overload


T = TypeVar("T")


class NamedRepository(Generic[T]):
    """Ordered repository addressed either by position or semantic name."""

    def __init__(self) -> None:
        self._items: list[T] = []
        self._names: dict[str, T] = {}

    def add(self, item: T) -> T:
        name = item.name
        if name in self._names:
            raise ValueError(f"duplicate name: {name}")
        self._items.append(item)
        self._names[name] = item
        return item

    def remove(self, key: int | str) -> T:
        item = self[key]
        self._items.remove(item)
        del self._names[item.name]
        return item

    def clear(self) -> None:
        self._items.clear()
        self._names.clear()

    @overload
    def __getitem__(self, key: int) -> T: ...

    @overload
    def __getitem__(self, key: str) -> T: ...

    def __getitem__(self, key: int | str) -> T:
        if isinstance(key, int):
            try:
                return self._items[key]
            except IndexError as exc:
                raise IndexError(f"repository index out of range: {key}") from exc
        if isinstance(key, str):
            try:
                return self._names[key]
            except KeyError as exc:
                raise KeyError(f"unknown name: {key}") from exc
        raise TypeError("repository key must be int or str")

    def __delitem__(self, key: int | str) -> None:
        self.remove(key)

    def index(self, key: str | T) -> int:
        item = self[key] if isinstance(key, str) else key
        return self._items.index(item)

    def names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self._items)

    def values(self) -> tuple[T, ...]:
        return tuple(self._items)

    def __contains__(self, key: object) -> bool:
        if isinstance(key, int):
            return -len(self._items) <= key < len(self._items)
        if isinstance(key, str):
            return key in self._names
        return False

    def __iter__(self) -> Iterator[T]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
