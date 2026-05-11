"""Eigenfrequency analysis step."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.supports import SupportCollector

from .constraint_method import ConstraintMethod


@dataclass(frozen=True, slots=True)
class ModalStep:
    name: str
    supports: tuple[SupportCollector, ...] = ()
    number_of_modes: int = 10
    constraint_method: ConstraintMethod | None = ConstraintMethod.NULLSPACE
