"""Repository for named load collectors."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .load_collector import LoadCollector, LoadEntry


class LoadCollectorRepository:
    """Repository for named groups of reusable loads."""

    def __init__(self) -> None:
        self._collectors: dict[str, LoadCollector] = {}

    def add(self, collector: LoadCollector) -> LoadCollector:
        if not isinstance(collector, LoadCollector):
            raise TypeError("collector must be a LoadCollector")
        self._collectors[collector.name] = collector
        return collector

    def get(self, name: str) -> LoadCollector:
        key = normalize_name(name)
        try:
            return self._collectors[key]
        except KeyError as exc:
            raise KeyError(f"unknown load collector: {key}") from exc

    def __getitem__(self, name: str) -> LoadCollector:
        return self.get(name)

    def has(self, value: str | LoadCollector) -> bool:
        if isinstance(value, LoadCollector):
            return any(item is value for item in self._collectors.values())
        return normalize_name(value) in self._collectors

    def all(self) -> tuple[LoadCollector, ...]:
        return tuple(self._collectors[key] for key in sorted(self._collectors))

    def __iter__(self) -> Iterator[LoadCollector]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._collectors)
