"""Surface definition data object."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SurfaceDefinition:
    surface_set: str
    target: object
    side: str | int
    id: int | None = None

    @property
    def explicit(self) -> bool:
        return self.id is not None
