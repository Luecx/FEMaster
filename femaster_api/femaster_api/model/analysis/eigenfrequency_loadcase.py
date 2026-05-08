"""Eigenfrequency loadcase definition."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod


@dataclass(frozen=True, slots=True)
class EigenfrequencyLoadcase:
    """Eigenfrequency extraction loadcase."""

    name: str
    supports: tuple[SupportCollector, ...] = ()
    number_of_modes: int = 10
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE

    def __post_init__(self) -> None:
        object.__setattr__(self, "supports", tuple(self.supports))
