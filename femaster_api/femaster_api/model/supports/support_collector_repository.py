"""Repository for named support collectors."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .support import Support
from .support_collector import SupportCollector


class SupportCollectorRepository:
    """Repository for named groups of reusable supports."""

    def __init__(self) -> None:
        self._collectors: dict[str, SupportCollector] = {}

    def add(self, collector: SupportCollector) -> SupportCollector:
        if not isinstance(collector, SupportCollector):
            raise TypeError("collector must be a SupportCollector")
        self._collectors[collector.name] = collector
        return collector

    def get(self, name: str) -> SupportCollector:
        key = normalize_name(name)
        try:
            return self._collectors[key]
        except KeyError as exc:
            raise KeyError(f"unknown support collector: {key}") from exc

    def __getitem__(self, name: str) -> SupportCollector:
        return self.get(name)

    def has(self, value: str | SupportCollector) -> bool:
        if isinstance(value, SupportCollector):
            return any(item is value for item in self._collectors.values())
        return normalize_name(value) in self._collectors

    def all(self) -> tuple[SupportCollector, ...]:
        return tuple(self._collectors[key] for key in sorted(self._collectors))

    def __iter__(self) -> Iterator[SupportCollector]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._collectors)
