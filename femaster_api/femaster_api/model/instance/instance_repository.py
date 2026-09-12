"""Repository of assembly instances."""

from ..common.named_repository import NamedRepository
from .instance import Instance


class InstanceRepository(NamedRepository[Instance]):
    """Named assembly-instance repository."""

    def export(self) -> str:
        return "\n\n".join(instance.export() for instance in self)
