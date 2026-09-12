"""Base Part-local section assignment."""

from ..common.named_object import NamedObject


class Section(NamedObject):
    """Base named section/property assignment."""

    def __init__(self, name: str, element_region: str) -> None:
        super().__init__(name)
        self.element_region = str(element_region)

    def export(self) -> str:
        raise NotImplementedError
