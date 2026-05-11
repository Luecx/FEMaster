"""Named reusable collection of supports."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .support import Support


@dataclass(frozen=True, slots=True)
class SupportCollector:
    name: str
    supports: tuple[Support, ...] = ()

    def add(self, support: Support) -> "SupportCollector":
        return replace(self, supports=(*self.supports, support))
