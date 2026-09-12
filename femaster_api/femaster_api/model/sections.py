"""Part-local FEMaster section assignments.

Sections are named Python objects for consistent repository access even when the
underlying FEMaster keyword does not require a NAME parameter. They reference
part-local element regions and global material/profile definitions by semantic
name, never by repository position.
"""

from __future__ import annotations

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository


class Section(NamedObject):
    """Base class for one named section assignment."""

    def __init__(self, name: str, element_region: str, material: str) -> None:
        super().__init__(name)
        self.element_region = str(element_region)
        self.material       = str(material)

    def to_femaster(self) -> str:
        raise NotImplementedError


class SolidSection(Section):
    """Three-dimensional continuum section assignment."""

    def __init__(self, name: str, element_region: str, material: str, orientation: str | None = None) -> None:
        super().__init__(name, element_region, material)
        self.orientation = orientation

    def to_femaster(self) -> str:
        return keyword(
            "SOLIDSECTION",
            ELSET=self.element_region,
            MATERIAL=self.material,
            ORIENTATION=self.orientation,
        )


class ShellSection(Section):
    """Homogeneous shell section with constant thickness."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        thickness: float,
        orientation: str | None = None,
    ) -> None:
        super().__init__(name, element_region, material)
        self.thickness   = float(thickness)
        self.orientation = orientation

    def to_femaster(self) -> str:
        return block([
            keyword(
                "SHELLSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                ORIENTATION=self.orientation,
            ),
            csv((self.thickness,)),
        ])


class BeamSection(Section):
    """Beam section referencing a named global profile."""

    def __init__(
        self,
        name: str,
        element_region: str,
        material: str,
        profile: str,
        orientation: tuple[float, float, float] | None = None,
    ) -> None:
        super().__init__(name, element_region, material)
        self.profile     = str(profile)
        self.orientation = orientation

    def to_femaster(self) -> str:
        lines = [
            keyword(
                "BEAMSECTION",
                ELSET=self.element_region,
                MATERIAL=self.material,
                PROFILE=self.profile,
            )
        ]
        if self.orientation is not None:
            lines.append(csv(self.orientation))
        return block(lines)


class TrussSection(Section):
    """Axial truss section defined by cross-sectional area."""

    def __init__(self, name: str, element_region: str, material: str, area: float) -> None:
        super().__init__(name, element_region, material)
        self.area = float(area)

    def to_femaster(self) -> str:
        return block([
            keyword("TRUSSSECTION", ELSET=self.element_region, MATERIAL=self.material),
            csv((self.area,)),
        ])


class SectionRepository(NamedRepository[Section]):
    """Named repository of part-local section assignments."""

    def to_femaster(self) -> str:
        return "\n\n".join(section.to_femaster() for section in self)
