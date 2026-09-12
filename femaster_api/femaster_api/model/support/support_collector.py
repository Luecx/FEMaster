"""Named owner of supports activated together by analysis steps.

Supports are stored directly in their collector rather than duplicated in a
global support repository.  This gives every support one clear ownership path
and mirrors the native loadcase syntax where steps activate collector names.
"""

from __future__ import annotations

from ..common.named_object import NamedObject
from .support import Support


class SupportCollector(NamedObject):
    """Named ordered collection of prescribed supports."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.supports: list[Support] = []

    def add(self, support: Support) -> Support:
        """Append one support and return the exact stored object."""

        if not isinstance(support, Support):
            raise TypeError("support must be a Support")
        self.supports.append(support)
        return support

    def export(self) -> str:
        """Export every support using this collector's semantic name."""

        return "\n\n".join(
            support.export(self.name)
            for support in self.supports
        )
