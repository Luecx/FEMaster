"""Collector-owned structural support definitions.

Supports prescribe nodal degrees of freedom and are owned directly by named
``SupportCollector`` objects.  Analysis steps activate collectors by semantic
name; the package therefore avoids a duplicate project-level repository of
individual support entries.
"""

from .support import Support
from .support_collector import SupportCollector
from .support_collector_repository import SupportCollectorRepository

__all__ = ["Support", "SupportCollector", "SupportCollectorRepository"]
