"""ABD shell elasticity data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ABDElasticity:
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.values) != 40:
            raise ValueError("ABD elasticity requires exactly 40 values")
