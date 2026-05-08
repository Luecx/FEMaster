"""Named containers for model object groups."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

from .entity_type import EntityType


@dataclass(frozen=True, slots=True)
class EntitySet:
    """Named group of model objects.

    Sets contain object references, not FEMaster IDs. Export IDs are resolved
    by the writer when a deck is serialized.
    """

    name: str
    members: tuple[object, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "members", tuple(self.members))

    @property
    def entity_type(self) -> EntityType:
        raise NotImplementedError

    def __getitem__(self, index: int | slice) -> object | tuple[object, ...]:
        if isinstance(index, slice):
            return self.members[index]
        return self.members[index]

    def __contains__(self, item: object) -> bool:
        return any(member is item for member in self.members)

    def __iter__(self) -> Iterator[object]:
        return iter(self.members)

    def __len__(self) -> int:
        return len(self.members)
