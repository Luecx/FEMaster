"""Repository of reusable parts with a permanent implicit default part.

FEMaster permits topology directly in root scope in addition to explicit
``*PART`` definitions.  The Python API represents that topology as a normal
``Part`` permanently stored at repository position zero.  There is deliberately
no second ``Project.default_part`` pointer or duplicated mesh state.

``default()`` therefore always returns ``self[0]`` and removal of that object is
forbidden.  Explicit user parts start at position one.
"""

from __future__ import annotations

from collections.abc import Iterator

from ..common.format import blocks
from ..common.named_repository import NamedRepository
from .part import Part


class PartRepository(NamedRepository[Part]):
    """Own all parts while protecting the implicit root/default part."""

    DEFAULT_NAME = "__DEFAULT_PART__"

    def __init__(self) -> None:
        super().__init__()
        super().add(Part(self.DEFAULT_NAME))

    def default(self) -> Part:
        """Return the permanent implicit part stored at index zero."""

        return self[0]

    def explicit(self) -> Iterator[Part]:
        """Iterate only explicitly named ``*PART`` definitions."""

        return iter(self._items[1:])

    def remove(self, key: int | str) -> Part:
        """Remove an explicit part while preventing default-part deletion."""

        part = self[key]
        if part is self.default():
            raise ValueError("the default part cannot be removed")
        return super().remove(key)

    def clear(self) -> None:
        """Remove all explicit parts and preserve the default part at index zero."""

        default = self.default()
        self._items = [default]
        self._names = {default.name: default}

    def export(self) -> str:
        """Export root topology first and explicit parts afterwards."""

        return blocks((
            self.default().export(root=True),
            *(part.export() for part in self.explicit()),
        ))
