from __future__ import annotations

from typing import Iterable, Iterator, List

from .constraint_base import ConstraintBase


class Constraints:
    """Container holding heterogeneous constraint entries."""

    def __init__(self, entries: Iterable[ConstraintBase] = ()) -> None:
        self._items: List[ConstraintBase] = []
        self.extend(entries)

    def add(self, entry: ConstraintBase) -> ConstraintBase:
        if not isinstance(entry, ConstraintBase):
            raise TypeError("Constraints can only contain ConstraintBase instances.")
        self._items.append(entry)
        return entry

    def extend(self, entries: Iterable[ConstraintBase]) -> None:
        for entry in entries:
            self.add(entry)

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[ConstraintBase]:
        return iter(self._items)

    def __getitem__(self, index: int) -> ConstraintBase:
        return self._items[index]

    def to_femaster(self) -> str:
        blocks: List[str] = []
        for entry in self._items:
            block = entry.to_femaster()
            if block:
                blocks.append(block)
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("Constraints.to_asami not implemented yet.")
