"""Repository for surface definitions."""

from __future__ import annotations

from typing import Iterator

from .surface_definition import SurfaceDefinition


class SurfaceRepository:
    """Repository for surface commands grouped by surface set."""

    def __init__(self) -> None:
        self._surfaces: list[SurfaceDefinition] = []

    def add(self, surface: SurfaceDefinition) -> SurfaceDefinition:
        if not isinstance(surface, SurfaceDefinition):
            raise TypeError("surface must be a SurfaceDefinition")
        self._surfaces.append(surface)
        return surface

    def all(self) -> tuple[SurfaceDefinition, ...]:
        return tuple(self._surfaces)

    def __contains__(self, surface: object) -> bool:
        return any(item is surface for item in self._surfaces)

    def sets(self) -> tuple[str, ...]:
        return tuple(sorted({surface.surface_set for surface in self._surfaces}))

    def by_set(self, surface_set: str) -> tuple[SurfaceDefinition, ...]:
        return tuple(surface for surface in self._surfaces if surface.surface_set == surface_set)

    def __iter__(self) -> Iterator[SurfaceDefinition]:
        return iter(self._surfaces)

    def __getitem__(self, index: int | slice) -> SurfaceDefinition | tuple[SurfaceDefinition, ...]:
        if isinstance(index, slice):
            return tuple(self._surfaces[index])
        return self._surfaces[index]

    def __len__(self) -> int:
        return len(self._surfaces)
