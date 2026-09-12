"""Named repository of assembly instances.

The repository emits only the individual ``*INSTANCE`` definitions.  ``Project``
owns the surrounding ``*ASSEMBLY`` scope because assembly-level regions,
surfaces and point-element properties must share that same block.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .instance import Instance


class InstanceRepository(NamedRepository[Instance]):
    """Own named assembly instances in deterministic creation order."""

    def export(self) -> str:
        """Export all instance blocks without adding an assembly wrapper."""

        return "\n\n".join(instance.export() for instance in self)
