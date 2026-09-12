"""Named repository of support collectors.

Steps reference support collectors through stable names.  Empty collectors are
retained in memory but omitted from native export until at least one support is
owned.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .support_collector import SupportCollector


class SupportCollectorRepository(NamedRepository[SupportCollector]):
    """Own named support collectors in deterministic order."""

    def export(self) -> str:
        return "\n\n".join(
            collector.export()
            for collector in self
            if collector.supports
        )
