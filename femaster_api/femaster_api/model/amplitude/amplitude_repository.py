"""Named repository of project-level amplitudes.

Amplitude names form the persistent references used by loads.  The repository
therefore provides deterministic name lookup and delegates all syntax to each
``Amplitude`` object.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .amplitude import Amplitude


class AmplitudeRepository(NamedRepository[Amplitude]):
    """Own globally named amplitude curves."""

    def export(self) -> str:
        return "\n\n".join(amplitude.export() for amplitude in self)
