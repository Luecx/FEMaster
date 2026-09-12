"""Restricted section repository for assembly-level point properties.

FEMaster currently permits concentrated mass, rotary inertia and spring
properties directly in assembly scope.  Continuum, shell, beam and truss section
assignments remain part-local.  This repository enforces that scope distinction
at insertion time rather than allowing an invalid project to survive until
export.
"""

from __future__ import annotations

from .section import Section
from .section_mass import MassSection
from .section_repository import SectionRepository
from .section_rotary_inertia import RotaryInertiaSection
from .section_spring import SpringSection


class AssemblySectionRepository(SectionRepository):
    """Assembly repository restricted to supported point-element properties."""

    _allowed = (MassSection, RotaryInertiaSection, SpringSection)

    def add(self, item: Section) -> Section:
        if not isinstance(item, self._allowed):
            raise TypeError(
                "assembly sections support only MassSection, "
                "RotaryInertiaSection and SpringSection"
            )
        return super().add(item)
