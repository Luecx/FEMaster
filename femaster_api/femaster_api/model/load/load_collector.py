"""Named load collector."""

from ..common.named_object import NamedObject
from .load import Load


class LoadCollector(NamedObject):
    """Named collection of loads activated together by analysis steps."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.loads: list[Load] = []

    def add(self, load: Load) -> Load:
        if not isinstance(load, Load):
            raise TypeError("load must derive from Load")
        self.loads.append(load)
        return load

    def export(self) -> str:
        return "\n\n".join(load.export(self.name) for load in self.loads)
