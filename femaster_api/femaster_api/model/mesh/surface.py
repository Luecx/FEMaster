"""Named element-boundary surface."""

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from ..common.typing import ElementReference


class Surface(NamedObject):
    """Named boundary surface valid in Part or Assembly scope."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.entries: list[tuple[ElementReference, int]] = []

    def add(self, element: ElementReference, side: int) -> "Surface":
        self.entries.append((element, int(side)))
        return self

    def export(self) -> str:
        return block([
            keyword("SURFACE", NAME=self.name, TYPE="ELEMENT"),
            *(csv((element, f"S{side}")) for element, side in self.entries),
        ])
