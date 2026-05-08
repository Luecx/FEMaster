"""Repository for named support collectors."""

from __future__ import annotations

from typing import Iterator

from femaster_api.model.names import require_name

from .support import Support
from .support_collector import SupportCollector


class SupportCollectorRepository:
    """Repository for named groups of reusable supports."""

    def __init__(self) -> None:
        self._collectors: dict[str, SupportCollector] = {}

    def add(self, collector: SupportCollector) -> SupportCollector:
        if not isinstance(collector, SupportCollector):
            raise TypeError("collector must be a SupportCollector")
        item = SupportCollector(require_name(collector.name), tuple(collector.supports))
        self._collectors[item.name] = item
        return item

    def get(self, name: str) -> SupportCollector:
        key = require_name(name)
        try:
            return self._collectors[key]
        except KeyError as exc:
            raise KeyError(f"unknown support collector: {key}") from exc

    def __getitem__(self, name: str) -> SupportCollector:
        return self.get(name)

    def has(self, value: str | SupportCollector) -> bool:
        if isinstance(value, SupportCollector):
            return any(item is value for item in self._collectors.values())
        return require_name(value) in self._collectors

    def __contains__(self, value: object) -> bool:
        if isinstance(value, str):
            return require_name(value) in self._collectors
        if isinstance(value, SupportCollector):
            return self.has(value)
        return False

    def all(self) -> tuple[SupportCollector, ...]:
        return tuple(self._collectors[key] for key in sorted(self._collectors))

    def __iter__(self) -> Iterator[SupportCollector]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._collectors)
