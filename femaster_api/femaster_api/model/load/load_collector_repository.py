"""Named repository of load collectors.

Analysis steps reference collectors by name, not by repository position.  Empty
collectors remain valid in memory but produce no deck output until they contain
at least one concrete load.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .load_collector import LoadCollector


class LoadCollectorRepository(NamedRepository[LoadCollector]):
    """Own named load collectors in deterministic order."""

    def export(self) -> str:
        return "\n\n".join(
            collector.export()
            for collector in self
            if collector.loads
        )
