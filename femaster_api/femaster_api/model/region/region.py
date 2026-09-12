"""Base named FEMaster region."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from ..common.typing import EntityReference


class Region(NamedObject):
    """Base region storing semantic entity references."""

    keyword_name: str | None = None

    def __init__(self, name: str, members: Iterable[EntityReference] = ()) -> None:
        super().__init__(name)
        self.members: list[EntityReference] = list(members)

    def add(self, *members: EntityReference) -> "Region":
        self.members.extend(members)
        return self

    def export(self) -> str:
        if self.keyword_name is None:
            raise NotImplementedError(
                f"{type(self).__name__} has no standalone native keyword"
            )

        rows = [
            csv(self.members[start:start + 16])
            for start in range(0, len(self.members), 16)
        ]
        return block([keyword(self.keyword_name, NAME=self.name), *rows])
