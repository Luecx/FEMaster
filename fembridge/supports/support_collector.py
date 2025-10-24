from __future__ import annotations

from typing import Iterable, List, Union

from .support_base import SupportBase
from .support import Support

SupportEntry = Union[SupportBase, Support]


class SupportCollector:
    """Collects supports and lets each write its FEMaster block."""

    def __init__(self, name: str, entries: Union[Iterable[SupportEntry], SupportEntry] = ()):
        if not name:
            raise ValueError("SupportCollector name must be non-empty.")
        self.name: str = str(name)

        if isinstance(entries, SupportBase):
            initial_entries = [entries]
        else:
            initial_entries = list(entries)

        self._items: List[SupportBase] = []
        self.extend(initial_entries)

    def add(self, support: SupportEntry) -> SupportEntry:
        if not isinstance(support, SupportBase):
            raise TypeError("SupportCollector can only contain SupportBase instances.")
        self._items.append(support)
        return support

    def extend(self, entries: Iterable[SupportEntry]) -> None:
        for support in entries:
            self.add(support)

    def __getitem__(self, idx: int) -> SupportBase:
        return self._items[idx]

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        blocks: List[str] = []
        for support in self._items:
            blocks.append(support.to_femaster(self.name))
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("SupportCollector.to_asami not implemented yet.")
