"""Part-local sections and assembly-level point-property assignments.

Every physical section/property form has its own ``section_*`` module.  Normal
``SectionRepository`` objects belong to parts, while
``AssemblySectionRepository`` deliberately permits only point-element
properties that FEMaster accepts directly in assembly scope.
"""

from .section import Section
from .section_beam import BeamSection
from .section_mass import MassSection
from .section_material import MaterialSection
from .section_repository import SectionRepository
from .section_repository_assembly import AssemblySectionRepository
from .section_rotary_inertia import RotaryInertiaSection
from .section_shell import ShellSection
from .section_shell_abd import ABDShellSection
from .section_solid import SolidSection
from .section_spring import SpringSection
from .section_truss import TrussSection

__all__ = [
    "ABDShellSection",
    "AssemblySectionRepository",
    "BeamSection",
    "MassSection",
    "MaterialSection",
    "RotaryInertiaSection",
    "Section",
    "SectionRepository",
    "ShellSection",
    "SolidSection",
    "SpringSection",
    "TrussSection",
]
