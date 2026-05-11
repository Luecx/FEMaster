"""Repository for reusable support definitions."""

from __future__ import annotations

from typing import Iterator

from .support import Support


class SupportRepository:
    """Repository for reusable support definitions."""

    def __init__(self) -> None:
        self._supports: list[Support] = []

    def add(self, support: Support) -> Support:
        if not isinstance(support, Support):
            raise TypeError("support must be a Support")
        self._supports.append(support)
        return support

    def all(self) -> tuple[Support, ...]:
        return tuple(self._supports)

    def __getitem__(self, index: int | slice) -> Support | tuple[Support, ...]:
        if isinstance(index, slice):
            return tuple(self._supports[index])
        return self._supports[index]

    def __iter__(self) -> Iterator[Support]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._supports)
