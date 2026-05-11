"""Base geometry entity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True, slots=True)
class GeometryEntity:
    tags: tuple[str, ...] = ()

    def contains(self, point: Iterable[float], *, tolerance: float = 1.0e-6) -> bool:
        raise NotImplementedError
