"""Repository of named support collectors."""

from ..common.named_repository import NamedRepository
from .support_collector import SupportCollector


class SupportCollectorRepository(NamedRepository[SupportCollector]):
    """Named support-collector repository."""

    def export(self) -> str:
        return "\n\n".join(
            collector.export() for collector in self if collector.supports
        )
