"""Part-local FEMaster section and point-property assignments.

Sections are named Python objects for consistent repository access even when the
underlying FEMaster keyword has no NAME parameter. They reference element regions
and shared definitions by semantic name, never by repository position.
"""

from __future__ import annotations

from typing import Iterable

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository


class Section(NamedObject):
    """Base class for one named section/property assignment."""

    def __init__(self, name: str, element_region: str) -> None:
        super().__init__(name)
        self.element_region = str(element_region)

    def export(self) -> str:
        raise NotImplementedError


class MaterialSection(Section):
    """Base class for sections referencing a global material definition."""

    def __init__(self, name: str, element_region: str, material: str) -> None:
        super().__init__(name, element_region)
        self.material = str(material)


class SolidSection(MaterialSection):
    """Three-dimensional continuum section assignment."""

    def __init__(self, name: str, element_region: str, material: str, orientation: str | None = None) -> None:
        super().__init__(name, element_region, material)
        self.orientation = orientation

    def export(self) -> str:
        return keyword(
            "SOLIDSECTION",
            ELSET=self.element_region,
            MATERIAL=self.material,
            ORIENTATION=self.orientation,
        )


class ShellSection(MaterialSection):
    """Material-integrated shell section with constant thickness."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        thickness: float,
        orientation: str | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region, material)
        self.thickness   = float(thickness)
        self.orientation = orientation
        self.csys_axis   = int(csys_axis)

        if self.csys_axis not in (1, 2, 3):
            raise ValueError("ShellSection csys_axis must be 1, 2 or 3")

    def export(self) -> str:
        return block([
            keyword(
                "SHELLSECTION",
                TYPE="INTEGRATED",
                ELSET=self.element_region,
                MATERIAL=self.material,
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            csv((self.thickness,)),
        ])


class ABDShellSection(Section):
    """Direct generalized shell stiffness using a 6x6 ABD and 2x2 shear matrix."""

    def __init__(
        self,
        name: str,
        element_region: str,
        abd: Iterable[float],
        shear: Iterable[float],
        *,
        thickness: float = 1.0,
        material: str | None = None,
        orientation: str | None = None,
        csys_axis: int = 1,
    ) -> None:
        super().__init__(name, element_region)
        self.abd         = tuple(float(value) for value in abd)
        self.shear       = tuple(float(value) for value in shear)
        self.thickness   = float(thickness)
        self.material    = material
        self.orientation = orientation
        self.csys_axis   = int(csys_axis)

        if len(self.abd) != 36:
            raise ValueError("ABDShellSection requires exactly 36 ABD values")
        if len(self.shear) != 4:
            raise ValueError("ABDShellSection requires exactly 4 shear values")
        if self.csys_axis not in (1, 2, 3):
            raise ValueError("ABDShellSection csys_axis must be 1, 2 or 3")

    def export(self) -> str:
        values = (*self.abd, *self.shear)
        rows = [csv(values[start:start + 8]) for start in range(0, len(values), 8)]
        return block([
            keyword(
                "SHELLSECTION",
                TYPE="ABD",
                ELSET=self.element_region,
                MATERIAL=self.material,
                THICKNESS=self.thickness,
                ORIENTATION=self.orientation,
                CSYSAXIS=self.csys_axis,
            ),
            *rows,
        ])


class BeamSection(MaterialSection):
    """Beam section referencing a named global profile and n1 direction."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        profile: str,
        orientation: Iterable[float],
    ) -> None:
        super().__init__(name, element_region, material)
        self.profile     = str(profile)
        self.orientation = tuple(float(value) for value in orientation)

        if len(self.orientation) != 3:
            raise ValueError("BeamSection orientation requires exactly 3 values")
        if all(value == 0.0 for value in self.orientation):
            raise ValueError("BeamSection orientation must be non-zero")

    def export(self) -> str:
        return block([
            keyword(
                "BEAMSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                PROFILE=self.profile,
            ),
            csv(self.orientation),
        ])


class TrussSection(MaterialSection):
    """Axial truss section defined by cross-sectional area."""

    def __init__(self, name: str, element_region: str, material: str, area: float) -> None:
        super().__init__(name, element_region, material)
        self.area = float(area)
        if self.area <= 0.0:
            raise ValueError("TrussSection area must be positive")

    def export(self) -> str:
        return block([
            keyword("TRUSSSECTION", ELSET=self.element_region, MATERIAL=self.material),
            csv((self.area,)),
        ])


class MassSection(Section):
    """Isotropic translational mass assigned to MASS point elements."""

    def __init__(self, name: str, element_region: str, mass: float) -> None:
        super().__init__(name, element_region)
        self.mass = float(mass)

    def export(self) -> str:
        return block([
            keyword("MASS", ELSET=self.element_region, TYPE="ISOTROPIC"),
            csv((self.mass,)),
        ])


class RotaryInertiaSection(Section):
    """Diagonal concentrated rotary inertia assigned to ROTARYI elements."""

    def __init__(self, name: str, element_region: str, inertia: Iterable[float]) -> None:
        super().__init__(name, element_region)
        self.inertia = tuple(float(value) for value in inertia)
        if len(self.inertia) != 3:
            raise ValueError("RotaryInertiaSection requires exactly 3 diagonal moments")

    def export(self) -> str:
        return block([
            keyword("ROTARY INERTIA", ELSET=self.element_region),
            csv((*self.inertia, 0.0, 0.0, 0.0)),
        ])


class SpringSection(Section):
    """One constant ground stiffness assigned to SPRING1 point elements."""

    def __init__(self, name: str, element_region: str, dof: int, stiffness: float) -> None:
        super().__init__(name, element_region)
        self.dof       = int(dof)
        self.stiffness = float(stiffness)
        if self.dof < 1 or self.dof > 6:
            raise ValueError("SpringSection dof must be between 1 and 6")

    def export(self) -> str:
        return block([
            keyword("SPRING", ELSET=self.element_region),
            csv((self.dof,)),
            csv((self.stiffness,)),
        ])


class SectionRepository(NamedRepository[Section]):
    """Named repository of Part-local section and point-property assignments."""

    def export(self) -> str:
        return "\n\n".join(section.export() for section in self)


class AssemblySectionRepository(SectionRepository):
    """Assembly-level point-property assignments valid after Part compilation.

    FEMaster currently permits MASS, ROTARY INERTIA and SPRING assignments in
    ASSEMBLY scope. Continuum, shell, beam and truss sections remain Part-local.
    """

    _allowed = (MassSection, RotaryInertiaSection, SpringSection)

    def add(self, item: Section) -> Section:
        if not isinstance(item, self._allowed):
            raise TypeError(
                "assembly sections currently support only MassSection, "
                "RotaryInertiaSection and SpringSection"
            )
        return super().add(item)
