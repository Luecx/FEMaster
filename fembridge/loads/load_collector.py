from __future__ import annotations

from typing import Iterable, List, Union

from .load_base import Load
from .pload import PLoad
from .dload import DLoad
from .cload import CLoad
from .vload import VLoad

LoadEntry = Union[PLoad, DLoad, CLoad, VLoad]


class LoadCollector:
    """Collects load entries and lets each entry produce its FEMaster block."""

    def __init__(self, name: str, entries: Union[Iterable[LoadEntry], LoadEntry] = ()):
        if not name:
            raise ValueError("LoadCollector name must be non-empty.")
        self.name: str = str(name)

        if isinstance(entries, Load):
            initial_entries = [entries]
        else:
            initial_entries = list(entries)

        self._items: List[Load] = []
        self.extend(initial_entries)

    def add(self, entry: LoadEntry) -> LoadEntry:
        if not isinstance(entry, Load):
            raise TypeError("LoadCollector can only contain Load instances.")
        self._items.append(entry)
        return entry

    def extend(self, entries: Iterable[LoadEntry]) -> None:
        for e in entries:
            self.add(e)

    def __getitem__(self, idx: int) -> Load:
        return self._items[idx]

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        blocks: List[str] = []
        for entry in self._items:
            blocks.append(entry.to_femaster(self.name))
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("LoadCollector.to_asami not implemented yet.")
