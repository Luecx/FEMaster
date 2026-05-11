"""Repository for section assignments."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .beam_section import BeamSection
from .shell_section import ShellSection
from .solid_section import SolidSection
from .truss_section import TrussSection

Vector3 = tuple[float, float, float]
Section = SolidSection | ShellSection | BeamSection | TrussSection


class SectionRepository:
    """Repository for section assignments."""

    def __init__(self) -> None:
        self._sections: dict[str, Section] = {}

    def add(self, section: Section) -> Section:
        if not isinstance(section, (SolidSection, ShellSection, BeamSection, TrussSection)):
            raise TypeError("section must be a section object")
        self._sections[normalize_name(section.name)] = section
        return section

    def get(self, name: str) -> Section:
        key = normalize_name(name)
        try:
            return self._sections[key]
        except KeyError as exc:
            raise KeyError(f"unknown section: {key}") from exc

    def __getitem__(self, name: str) -> Section:
        return self.get(name)

    def has(self, name: str) -> bool:
        return normalize_name(name) in self._sections

    def all(self) -> tuple[Section, ...]:
        return tuple(self._sections[key] for key in sorted(self._sections))

    def __iter__(self) -> Iterator[Section]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._sections)
