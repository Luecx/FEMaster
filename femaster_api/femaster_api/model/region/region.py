"""Generic base for named regions that own references to real model objects.

A region groups already-existing entities and therefore stores those entities
as Python objects rather than integer IDs or semantic-name strings.  Concrete
region classes specialize the member type and define how their objects are
reduced to native IDs/names during export.  The generic base only owns ordering,
validation and common keyword formatting.

This distinction is important for the public API: a ``NodeRegion`` contains
``Node`` instances, an ``ElementRegion`` contains ``Element`` instances, and a
``SurfaceRegion`` contains ``Surface`` instances.  Parsing is responsible for
resolving native tokens to those objects before a region is constructed.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Generic, TypeVar

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject


TRegionMember = TypeVar("TRegionMember")


class Region(NamedObject, Generic[TRegionMember]):
    """Base class for one ordered named region of concrete model objects."""

    keyword_name: str | None = None
    name_key = "NAME"
    member_type: type[object] = object

    def __init__(
        self,
        name: str,
        members: Iterable[TRegionMember] = (),
    ) -> None:
        super().__init__(name)
        self.members: list[TRegionMember] = []
        self.add(*members)

    def add(self, *members: TRegionMember) -> "Region[TRegionMember]":
        """Append validated object members and return this region."""

        for member in members:
            if not isinstance(member, self.member_type):
                raise TypeError(
                    f"{type(self).__name__} members must be "
                    f"{self.member_type.__name__} objects"
                )
        self.members.extend(members)
        return self

    def _member_value(self, member: TRegionMember) -> int | str:
        """Return the native ID/name used to serialize one member."""

        raise NotImplementedError

    def export(self) -> str:
        """Export a region that has a direct native set keyword."""

        if self.keyword_name is None:
            raise NotImplementedError(
                f"{type(self).__name__} has no standalone FEMaster keyword"
            )

        values = [self._member_value(member) for member in self.members]
        rows = [
            csv(values[start:start + 16])
            for start in range(0, len(values), 16)
        ]
        return block([
            keyword(self.keyword_name, **{self.name_key: self.name}),
            *rows,
        ])
