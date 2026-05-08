"""Node objects used by the in-memory model."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True, eq=False)
class Node:
    """A geometric finite-element node in global coordinates.

    Nodes intentionally do not store FEMaster export IDs. The writer assigns
    local IDs only while serializing a deck.
    """

    x: float
    y: float
    z: float
