"""Tie constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import SurfaceSet


@dataclass(frozen=True, slots=True)
class TieConstraint:
    """Tie two surface sets together."""

    master: SurfaceSet
    slave: SurfaceSet
    distance: float
    adjust: bool = False
