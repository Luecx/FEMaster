"""Repository for object-backed entity sets."""

from __future__ import annotations

from dataclasses import replace
from typing import Iterable, Iterator

from femaster_api.model.elements import Element
from femaster_api.model.names import require_name
from femaster_api.model.nodes import Node

from .element_set import ElementSet
from .entity_set import EntitySet
from .entity_type import EntityType
from .member_ids import unique_members
from .node_set import NodeSet
from .surface_set import SurfaceSet


class SetRepository:
    """Repository for all FEM entity sets."""

    def __init__(self) -> None:
        self._sets: dict[EntityType, dict[str, EntitySet]] = {
            EntityType.NODE: {},
            EntityType.ELEMENT: {},
            EntityType.SURFACE: {},
        }

    def add(self, entity_set: EntitySet) -> EntitySet:
        if not isinstance(entity_set, EntitySet):
            raise TypeError("entity_set must be an EntitySet")
        _validate_members(entity_set)
        item = replace(
            entity_set,
            name=require_name(entity_set.name),
            members=unique_members(entity_set.members),
        )
        self._sets[item.entity_type][item.name] = item
        return item

    def get(self, entity_type: EntityType, name: str) -> EntitySet:
        set_name = require_name(name)
        try:
            return self._sets[entity_type][set_name]
        except KeyError as exc:
            raise KeyError(f"unknown {entity_type.value} set: {set_name}") from exc

    def __getitem__(self, key: tuple[EntityType, str] | tuple[str, str]) -> EntitySet:
        entity_type, name = key
        if not isinstance(entity_type, EntityType):
            entity_type = EntityType(str(entity_type).lower())
        return self.get(entity_type, name)

    def has(self, value: str | EntitySet, entity_type: EntityType | None = None) -> bool:
        if isinstance(value, EntitySet):
            if entity_type is not None and value.entity_type != entity_type:
                return False
            return any(item is value for item in self._sets[value.entity_type].values())
        if entity_type is None:
            raise ValueError("entity_type is required when checking a set by name")
        set_name = require_name(value)
        return set_name in self._sets[entity_type]

    def all(self, entity_type: EntityType | None = None) -> tuple[EntitySet, ...]:
        if entity_type is not None:
            values: Iterable[EntitySet] = self._sets[entity_type].values()
        else:
            values = (
                item
                for entity_sets in self._sets.values()
                for item in entity_sets.values()
            )
        return tuple(sorted(values, key=lambda item: (item.entity_type.value, item.name)))

    def __iter__(self) -> Iterator[EntitySet]:
        return iter(self.all())

    def __len__(self) -> int:
        return sum(len(items) for items in self._sets.values())


def _validate_members(entity_set: EntitySet) -> None:
    for member in entity_set.members:
        if isinstance(entity_set, NodeSet) and not isinstance(member, Node):
            raise TypeError("NodeSet members must be Node objects")
        if isinstance(entity_set, ElementSet) and not isinstance(member, Element):
            raise TypeError("ElementSet members must be Element objects")
