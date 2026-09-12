"""Repository of named load collectors."""

from ..common.named_repository import NamedRepository
from .load_collector import LoadCollector


class LoadCollectorRepository(NamedRepository[LoadCollector]):
    """Named load-collector repository."""

    def export(self) -> str:
        return "\n\n".join(
            collector.export() for collector in self if collector.loads
        )
