"""Repository of global model fields."""

from ..common.named_repository import NamedRepository
from .field import Field


class FieldRepository(NamedRepository[Field]):
    """Named global field repository."""

    def export(self) -> str:
        return "\n\n".join(field.export() for field in self)
