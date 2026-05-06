"""Element set data object."""

from __future__ import annotations

from dataclasses import dataclass

from .entity_set import EntitySet
from .entity_type import EntityType


@dataclass(frozen=True, slots=True)
class ElementSet(EntitySet):
    @property
    def entity_type(self) -> EntityType:
        return EntityType.ELEMENT
