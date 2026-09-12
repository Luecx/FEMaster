"""Assembly-level point-property repository."""

from .mass_section import MassSection
from .rotary_inertia_section import RotaryInertiaSection
from .section import Section
from .section_repository import SectionRepository
from .spring_section import SpringSection


class AssemblySectionRepository(SectionRepository):
    """Assembly repository restricted to point-element properties."""

    def add(self, item: Section) -> Section:
        if not isinstance(item, (MassSection, RotaryInertiaSection, SpringSection)):
            raise TypeError(
                "assembly sections support only MassSection, "
                "RotaryInertiaSection and SpringSection"
            )
        return super().add(item)
