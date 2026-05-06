"""Node data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True, eq=False)
class Node:
    """A geometric FEM node in global coordinates."""

    x: float
    y: float
    z: float
    id: int | None = None
