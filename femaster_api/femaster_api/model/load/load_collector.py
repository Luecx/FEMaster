"""Named owner of loads activated together by analysis steps.

A collector owns its load objects directly.  This is intentionally stronger
than storing only references into a second global load repository: one load has
one clear lifetime and one collector context.  During export the collector name
is passed to each concrete load so the load can emit its own keyword.
"""

from __future__ import annotations

from ..common.named_object import NamedObject
from .load import Load


class LoadCollector(NamedObject):
    """Named ordered collection of concrete load objects."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.loads: list[Load] = []

    def add(self, load: Load) -> Load:
        """Append a load and return the exact stored object."""

        if not isinstance(load, Load):
            raise TypeError("load must derive from Load")
        self.loads.append(load)
        return load

    def export(self) -> str:
        """Export every owned load using this collector's semantic name."""

        return "\n\n".join(
            load.export(self.name)
            for load in self.loads
        )
