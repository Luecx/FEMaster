"""Surface set data object."""

from __future__ import annotations

from dataclasses import dataclass

from .entity_set import EntitySet
from .entity_type import EntityType


@dataclass(frozen=True, slots=True)
class SurfaceSet(EntitySet):
    """Named group of surface definitions."""

    @property
    def entity_type(self) -> EntityType:
        return EntityType.SURFACE
