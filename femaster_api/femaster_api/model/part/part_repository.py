"""Repository of reusable Parts with permanent default Part."""

from collections.abc import Iterator

from ..common.format import blocks
from ..common.named_repository import NamedRepository
from .part import Part


class PartRepository(NamedRepository[Part]):
    """Part repository whose index zero is permanently the implicit default Part."""

    DEFAULT_NAME = "__DEFAULT_PART__"

    def __init__(self) -> None:
        super().__init__()
        super().add(Part(self.DEFAULT_NAME))

    def default(self) -> Part:
        return self[0]

    def explicit(self) -> Iterator[Part]:
        return iter(self._items[1:])

    def remove(self, key: int | str) -> Part:
        part = self[key]
        if part is self.default():
            raise ValueError("the default part cannot be removed")
        return super().remove(key)

    def clear(self) -> None:
        default = self.default()
        self._items = [default]
        self._names = {default.name: default}

    def export(self) -> str:
        return blocks((
            self.default().export(root=True),
            *(part.export() for part in self.explicit()),
        ))
