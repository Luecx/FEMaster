"""Named repositories for object-backed entity sets."""

from __future__ import annotations

from dataclasses import replace
from typing import Generic, Iterable, Iterator, TypeVar

from femaster_api.model.elements import Element
from femaster_api.model.names import require_name
from femaster_api.model.nodes import Node
from femaster_api.model.surfaces import SurfaceDefinition

from .element_set import ElementSet
from .entity_set import EntitySet
from .entity_type import EntityType
from .member_ids import unique_members
from .node_set import NodeSet
from .surface_set import SurfaceSet

SetT = TypeVar("SetT", bound=EntitySet)


class NamedSetRepository(Generic[SetT]):
    """Named container for one kind of entity set."""

    entity_type: EntityType
    set_type: type[SetT]

    def __init__(self) -> None:
        self._items: dict[str, SetT] = {}

    def add(self, entity_set: SetT) -> SetT:
        """Add a set and return the stored object."""

        if not isinstance(entity_set, self.set_type):
            raise TypeError(f"entity_set must be a {self.set_type.__name__}")
        _validate_members(entity_set)
        item = replace(
            entity_set,
            name=require_name(entity_set.name),
            members=unique_members(entity_set.members),
        )
        self._items[item.name] = item
        return item

    def get(self, name: str) -> SetT:
        key = require_name(name)
        try:
            return self._items[key]
        except KeyError as exc:
            raise KeyError(f"unknown {self.entity_type.value} set: {key}") from exc

    def has(self, value: str | SetT) -> bool:
        if isinstance(value, self.set_type):
            return any(item is value for item in self._items.values())
        return require_name(value) in self._items

    def all(self) -> tuple[SetT, ...]:
        return tuple(self._items[key] for key in sorted(self._items))

    def __getitem__(self, name: str) -> SetT:
        return self.get(name)

    def __contains__(self, value: object) -> bool:
        if isinstance(value, str):
            return require_name(value) in self._items
        if isinstance(value, self.set_type):
            return self.has(value)
        return False

    def __iter__(self) -> Iterator[SetT]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._items)


class NodeSetRepository(NamedSetRepository[NodeSet]):
    """Named container for node sets."""

    entity_type = EntityType.NODE
    set_type = NodeSet


class ElementSetRepository(NamedSetRepository[ElementSet]):
    """Named container for element sets."""

    entity_type = EntityType.ELEMENT
    set_type = ElementSet


class SurfaceSetRepository(NamedSetRepository[SurfaceSet]):
    """Named container for surface sets."""

    entity_type = EntityType.SURFACE
    set_type = SurfaceSet


class SetRepository:
    """Compatibility repository for all FEM entity sets.

    New code should prefer ``model.node_sets``, ``model.element_sets``, and
    ``model.surface_sets``.
    """

    def __init__(self) -> None:
        self._repositories: dict[EntityType, NamedSetRepository] = {
            EntityType.NODE: NodeSetRepository(),
            EntityType.ELEMENT: ElementSetRepository(),
            EntityType.SURFACE: SurfaceSetRepository(),
        }

    def add(self, entity_set: EntitySet) -> EntitySet:
        if not isinstance(entity_set, EntitySet):
            raise TypeError("entity_set must be an EntitySet")
        return self._repositories[entity_set.entity_type].add(entity_set)

    def get(self, entity_type: EntityType, name: str) -> EntitySet:
        return self._repositories[entity_type].get(name)

    def __getitem__(self, key: tuple[EntityType, str] | tuple[str, str]) -> EntitySet:
        entity_type, name = key
        if not isinstance(entity_type, EntityType):
            entity_type = EntityType(str(entity_type).lower())
        return self.get(entity_type, name)

    def has(self, value: str | EntitySet, entity_type: EntityType | None = None) -> bool:
        if isinstance(value, EntitySet):
            if entity_type is not None and value.entity_type != entity_type:
                return False
            return self._repositories[value.entity_type].has(value)
        if entity_type is None:
            raise ValueError("entity_type is required when checking a set by name")
        return self._repositories[entity_type].has(value)

    def all(self, entity_type: EntityType | None = None) -> tuple[EntitySet, ...]:
        if entity_type is not None:
            values: Iterable[EntitySet] = self._repositories[entity_type]
        else:
            values = (
                item
                for repository in self._repositories.values()
                for item in repository
            )
        return tuple(sorted(values, key=lambda item: (item.entity_type.value, item.name)))

    def __contains__(self, value: object) -> bool:
        if isinstance(value, EntitySet):
            return self.has(value)
        return False

    def __iter__(self) -> Iterator[EntitySet]:
        return iter(self.all())

    def __len__(self) -> int:
        return sum(len(repository) for repository in self._repositories.values())


def _validate_members(entity_set: EntitySet) -> None:
    for member in entity_set.members:
        if isinstance(entity_set, NodeSet) and not isinstance(member, Node):
            raise TypeError("NodeSet members must be Node objects")
        if isinstance(entity_set, ElementSet) and not isinstance(member, Element):
            raise TypeError("ElementSet members must be Element objects")
        if isinstance(entity_set, SurfaceSet) and not isinstance(member, SurfaceDefinition):
            raise TypeError("SurfaceSet members must be SurfaceDefinition objects")
