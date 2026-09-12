"""Ordered repository for model objects with semantic names.

``NamedRepository`` intentionally supports two distinct access modes:
``repository[index]`` is insertion-order convenience access, whereas
``repository[name]`` is semantic lookup.  Positional indices are never written
into a FEMaster deck and must not be used as persistent cross-object references.

The generic type remains intentionally unconstrained at the type-variable level
so this module does not need a second protocol class.  Concrete repository users
are expected to provide objects exposing the ``name`` property established by
``NamedObject``.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Generic, TypeVar, overload


TNamed = TypeVar("TNamed")


class NamedRepository(Generic[TNamed]):
    """Own named objects in deterministic insertion order."""

    def __init__(self) -> None:
        self._items: list[TNamed] = []
        self._names: dict[str, TNamed] = {}

    def add(self, item: TNamed) -> TNamed:
        """Insert ``item`` and return the exact stored object."""

        name = item.name
        if name in self._names:
            raise ValueError(f"duplicate name: {name}")
        self._items.append(item)
        self._names[name] = item
        return item

    def remove(self, key: int | str) -> TNamed:
        """Remove and return an item selected by position or semantic name."""

        item = self[key]
        self._items.remove(item)
        del self._names[item.name]
        return item

    def clear(self) -> None:
        """Remove every owned object."""

        self._items.clear()
        self._names.clear()

    @overload
    def __getitem__(self, key: int) -> TNamed: ...

    @overload
    def __getitem__(self, key: str) -> TNamed: ...

    def __getitem__(self, key: int | str) -> TNamed:
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

    def index(self, key: str | TNamed) -> int:
        """Return the current insertion-order position of an item."""

        item = self[key] if isinstance(key, str) else key
        return self._items.index(item)

    def names(self) -> tuple[str, ...]:
        """Return all semantic names in deterministic order."""

        return tuple(item.name for item in self._items)

    def values(self) -> tuple[TNamed, ...]:
        """Return an immutable snapshot of owned objects."""

        return tuple(self._items)

    def __contains__(self, key: object) -> bool:
        if isinstance(key, int):
            return -len(self._items) <= key < len(self._items)
        if isinstance(key, str):
            return key in self._names
        return False

    def __iter__(self) -> Iterator[TNamed]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
