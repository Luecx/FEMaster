"""Named entity set data object."""

from __future__ import annotations

from dataclasses import dataclass

from .entity_type import EntityType


@dataclass(frozen=True, slots=True)
class EntitySet:
    name: str
    members: tuple[object, ...]
    generated: bool = False

    @property
    def entity_type(self) -> EntityType:
        raise NotImplementedError
