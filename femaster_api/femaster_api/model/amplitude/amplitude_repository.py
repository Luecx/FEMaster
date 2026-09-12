"""Repository of global amplitudes."""

from ..common.named_repository import NamedRepository
from .amplitude import Amplitude


class AmplitudeRepository(NamedRepository[Amplitude]):
    """Named global amplitude repository."""

    def export(self) -> str:
        return "\n\n".join(amplitude.export() for amplitude in self)
