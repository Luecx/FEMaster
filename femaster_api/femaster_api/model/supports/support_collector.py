"""Named reusable collection of supports."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .support import Support


@dataclass(frozen=True, slots=True)
class SupportCollector:
    """Named list of supports that can be attached to loadcases."""

    name: str
    supports: tuple[Support, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "supports", tuple(self.supports))

    def add(self, support: Support) -> "SupportCollector":
        return replace(self, supports=(*self.supports, support))
