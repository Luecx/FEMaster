"""Named support collector."""

from ..common.named_object import NamedObject
from .support import Support


class SupportCollector(NamedObject):
    """Named collection of prescribed supports."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.supports: list[Support] = []

    def add(self, support: Support) -> Support:
        if not isinstance(support, Support):
            raise TypeError("support must be a Support")
        self.supports.append(support)
        return support

    def export(self) -> str:
        return "\n\n".join(
            support.export(self.name) for support in self.supports
        )
