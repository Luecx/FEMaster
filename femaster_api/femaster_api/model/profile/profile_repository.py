"""Named repository of global beam profiles.

Profiles are shared project-level definitions referenced from part-local beam
sections by immutable semantic name.  Repository positions are convenience
access only and are never exported into the deck.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .profile import Profile


class ProfileRepository(NamedRepository[Profile]):
    """Own globally named beam profiles."""

    def export(self) -> str:
        return "\n\n".join(profile.export() for profile in self)
