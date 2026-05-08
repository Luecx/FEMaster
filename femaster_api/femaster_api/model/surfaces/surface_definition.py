"""Surface definition model object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SurfaceDefinition:
    """One element-side surface definition.

    Surface definitions do not store export IDs. Surface sets and the writer
    reference these objects directly.
    """

    surface_set: str
    target: object
    side: str | int


Surface = SurfaceDefinition
