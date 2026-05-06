"""Transient time control data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TimeControl:
    start: float
    end: float
    step: float
