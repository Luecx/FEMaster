"""Ordered repository primitives used by the FEMaster Python API.

The repositories in this module deliberately implement only ownership, lookup
and deterministic iteration. FEM semantics remain in the concrete model
objects. Named repositories support both positional and name-based access,
while ID repositories preserve sparse FEMaster identifiers exactly as supplied.
"""

from __future__ import annotations

from typing import Generic, Iterator, Protocol, TypeVar, overload


class Named(Protocol):
    """Protocol implemented by objects with a stable semantic name."""

    @property
    def name(self) -> str: ...


class Identified(Protocol):
    """Protocol implemented by objects carrying a FEMaster integer identifier."""

    id: int


TNamed = TypeVar("TNamed", bound=Named)
TId = TypeVar("TId", bound=Identified)


class NamedObject:
    """Base class for FEMaster objects with an immutable semantic name.

    Names are immutable because NamedRepository maintains a direct lookup table
    from name to object. Mutable names would invalidate that table silently.
    """

    __slots__ = ("_name",)

    def __init__(self, name: str) -> None:
        name = str(name).strip()
        if not name:
            raise ValueError("name must not be empty")
        self._name = name

    @property
    def name(self) -> str:
        """Return the immutable semantic name."""

        return self._name

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r})"


class NamedRepository(Generic[TNamed]):
    """Ordered repository for uniquely named objects.

    Integer keys address the current repository position. String keys address
    the persistent semantic name. Positions are convenience indices only and
    may change after removal; model references must therefore use names.
    """

    def __init__(self) -> None:
        self._items: list[TNamed] = []
        self._names: dict[str, TNamed] = {}

    # ------------------------------------------------------------------
    # Modification
    # ------------------------------------------------------------------

    def add(self, item: TNamed) -> TNamed:
        """Add an object while preserving insertion order."""

        if item.name in self._names:
            raise ValueError(f"duplicate name: {item.name}")

        self._items.append(item)
        self._names[item.name] = item
        return item

    def remove(self, key: int | str) -> TNamed:
        """Remove and return an object by position or semantic name."""

        item = self[key]
        self._items.remove(item)
        del self._names[item.name]
        return item

    def clear(self) -> None:
        """Remove every object from the repository."""

        self._items.clear()
        self._names.clear()

    # ------------------------------------------------------------------
    # Access
    # ------------------------------------------------------------------

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
        """Return the current positional index of an object."""

        item = self[key] if isinstance(key, str) else key
        return self._items.index(item)

    def names(self) -> tuple[str, ...]:
        """Return semantic names in deterministic repository order."""

        return tuple(item.name for item in self._items)

    def values(self) -> tuple[TNamed, ...]:
        """Return all objects in deterministic repository order."""

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


class IdRepository(Generic[TId]):
    """Ordered repository preserving sparse FEMaster integer identifiers.

    Integer subscription addresses the semantic FEMaster ID, not the insertion
    position. Positional access is intentionally explicit through at() so IDs
    and Python sequence positions cannot be confused.
    """

    def __init__(self) -> None:
        self._items: list[TId] = []
        self._ids: dict[int, TId] = {}

    def add(self, item: TId) -> TId:
        """Add an object without renumbering its FEMaster identifier."""

        if item.id in self._ids:
            raise ValueError(f"duplicate id: {item.id}")

        self._items.append(item)
        self._ids[item.id] = item
        return item

    def remove(self, id: int) -> TId:
        """Remove and return an object by FEMaster identifier."""

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
        """Return an object by current insertion-order position."""

        try:
            return self._items[index]
        except IndexError as exc:
            raise IndexError(f"repository index out of range: {index}") from exc

    def ids(self) -> tuple[int, ...]:
        """Return FEMaster identifiers in deterministic repository order."""

        return tuple(item.id for item in self._items)

    def values(self) -> tuple[TId, ...]:
        """Return all objects in deterministic repository order."""

        return tuple(self._items)

    def __contains__(self, id: object) -> bool:
        return isinstance(id, int) and id in self._ids

    def __iter__(self) -> Iterator[TId]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
