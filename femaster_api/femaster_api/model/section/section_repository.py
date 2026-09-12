"""Repository of Part-local sections."""

from ..common.named_repository import NamedRepository
from .section import Section


class SectionRepository(NamedRepository[Section]):
    """Named Part-local section/property repository."""

    def export(self) -> str:
        return "\n\n".join(section.export() for section in self)
