"""Newmark integration control data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class NewmarkControl:
    beta: float = 0.25
    gamma: float = 0.5
