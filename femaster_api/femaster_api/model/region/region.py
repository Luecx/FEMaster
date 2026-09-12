"""Base class for named FEMaster regions.

A region is a semantic set of existing entities; it does not create topology.
Concrete subclasses select the entity domain and, where FEMaster has a direct
keyword representation, the keyword used during export.  Members preserve their
input order and may be sparse integer IDs or semantic string references.

Topology definitions such as ``ElementSurface`` deliberately live in the
separate ``surface`` package.  ``SurfaceRegion`` only groups surfaces.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from ..common.typing import EntityReference


class Region(NamedObject):
    """Base class for one ordered named entity region."""

    keyword_name: str | None = None
    name_key = "NAME"

    def __init__(
        self,
        name: str,
        members: Iterable[EntityReference] = (),
    ) -> None:
        super().__init__(name)
        self.members: list[EntityReference] = list(members)

    def add(self, *members: EntityReference) -> "Region":
        """Append semantic member references and return this region."""

        self.members.extend(members)
        return self

    def export(self) -> str:
        """Export a region that has a direct native set keyword."""

        if self.keyword_name is None:
            raise NotImplementedError(
                f"{type(self).__name__} has no standalone FEMaster keyword"
            )

        rows = [
            csv(self.members[start:start + 16])
            for start in range(0, len(self.members), 16)
        ]
        return block([
            keyword(self.keyword_name, **{self.name_key: self.name}),
            *rows,
        ])
