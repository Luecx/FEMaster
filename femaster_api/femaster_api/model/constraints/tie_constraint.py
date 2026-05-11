"""Tie constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import EntitySet


@dataclass(frozen=True, slots=True)
class TieConstraint:
    master: EntitySet
    slave: EntitySet
    distance: float
    adjust: bool = False
