"""Repository for section assignments and beam profiles."""

from __future__ import annotations

from typing import Iterator

from femaster_api.model.names import require_name

from .beam_section import BeamSection
from .profile import Profile
from .shell_section import ShellSection
from .solid_section import SolidSection
from .truss_section import TrussSection

Vector3 = tuple[float, float, float]
Section = SolidSection | ShellSection | BeamSection | TrussSection


class SectionRepository:
    """Repository for profile and section assignments."""

    def __init__(self) -> None:
        self._sections: dict[str, Section] = {}
        self._profiles: dict[str, Profile] = {}

    def add(self, section: Section) -> Section:
        if not isinstance(section, (SolidSection, ShellSection, BeamSection, TrussSection)):
            raise TypeError("section must be a section object")
        self._sections[require_name(section.name)] = section
        return section

    def add_profile(self, profile: Profile) -> Profile:
        if not isinstance(profile, Profile):
            raise TypeError("profile must be a Profile")
        self._profiles[require_name(profile.name)] = profile
        return profile

    def get(self, name: str) -> Section:
        key = require_name(name)
        try:
            return self._sections[key]
        except KeyError as exc:
            raise KeyError(f"unknown section: {key}") from exc

    def __getitem__(self, name: str) -> Section:
        return self.get(name)

    def has(self, value: str | Section) -> bool:
        if isinstance(value, (SolidSection, ShellSection, BeamSection, TrussSection)):
            return any(item is value for item in self._sections.values())
        return require_name(value) in self._sections

    def __contains__(self, value: object) -> bool:
        if isinstance(value, str):
            return require_name(value) in self._sections
        if isinstance(value, (SolidSection, ShellSection, BeamSection, TrussSection)):
            return self.has(value)
        return False

    def profiles(self) -> tuple[Profile, ...]:
        return tuple(self._profiles[key] for key in sorted(self._profiles))

    def has_profile(self, value: str | Profile) -> bool:
        if isinstance(value, Profile):
            return any(item is value for item in self._profiles.values())
        return require_name(value) in self._profiles

    def all(self) -> tuple[Section, ...]:
        return tuple(self._sections[key] for key in sorted(self._sections))

    def __iter__(self) -> Iterator[Section]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._sections)
