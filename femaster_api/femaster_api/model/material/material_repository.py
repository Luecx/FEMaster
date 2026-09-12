"""Repository of global materials."""

from ..common.named_repository import NamedRepository
from .material import Material


class MaterialRepository(NamedRepository[Material]):
    """Named global material repository."""

    def export(self) -> str:
        return "\n\n".join(material.export() for material in self)
