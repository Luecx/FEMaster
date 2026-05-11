"""Section model objects."""

from .beam_section import BeamSection
from .section_repository import SectionRepository
from .shell_section import ShellSection
from .solid_section import SolidSection
from .truss_section import TrussSection

Section = SolidSection | ShellSection | BeamSection | TrussSection

__all__ = [
    "BeamSection",
    "Section",
    "SectionRepository",
    "ShellSection",
    "SolidSection",
    "TrussSection",
]
