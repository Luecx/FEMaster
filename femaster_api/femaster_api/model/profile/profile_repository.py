"""Repository of global beam profiles."""

from ..common.named_repository import NamedRepository
from .profile import Profile


class ProfileRepository(NamedRepository[Profile]):
    """Named global beam-profile repository."""

    def export(self) -> str:
        return "\n\n".join(profile.export() for profile in self)
