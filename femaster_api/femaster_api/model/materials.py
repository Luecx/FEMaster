"""Global FEMaster material and profile definitions.

Material behavior is represented by small explicit classes. Each constitutive
sub-definition knows how to export only its own keyword block; Material composes
those blocks in FEMaster order. No serializer registry is involved.
"""

from __future__ import annotations

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository


class Elasticity:
    """Base class for exportable elastic constitutive definitions."""

    def export(self) -> str:
        raise NotImplementedError


class IsotropicElasticity(Elasticity):
    """Linear isotropic elasticity defined by Young's modulus and Poisson ratio."""

    def __init__(self, youngs_modulus: float, poisson_ratio: float) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio  = float(poisson_ratio)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="ISOTROPIC"),
            csv((self.youngs_modulus, self.poisson_ratio)),
        ])


class GeneralizedIsotropicElasticity(Elasticity):
    """Generalized isotropic elasticity with an independent shear modulus."""

    def __init__(self, youngs_modulus: float, poisson_ratio: float, shear_modulus: float) -> None:
        self.youngs_modulus = float(youngs_modulus)
        self.poisson_ratio  = float(poisson_ratio)
        self.shear_modulus  = float(shear_modulus)

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="GENISO"),
            csv((self.youngs_modulus, self.poisson_ratio, self.shear_modulus)),
        ])


class OrthotropicElasticity(Elasticity):
    """Orthotropic engineering constants in material coordinates."""

    def __init__(
        self,
        e1: float,
        e2: float,
        e3: float,
        nu12: float,
        nu13: float,
        nu23: float,
        g12: float,
        g13: float,
        g23: float,
    ) -> None:
        self.values = tuple(float(value) for value in (e1, e2, e3, nu12, nu13, nu23, g12, g13, g23))

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="ENGINEERINGCONSTANTS"),
            csv(self.values),
        ])


class ABDElasticity(Elasticity):
    """General shell ABD constitutive matrix stored in FEMaster input order."""

    def __init__(self, values: tuple[float, ...] | list[float]) -> None:
        self.values = tuple(float(value) for value in values)

    def export(self) -> str:
        lines = [keyword("ELASTIC", TYPE="ABD")]
        for start in range(0, len(self.values), 8):
            lines.append(csv(self.values[start:start + 8]))
        return block(lines)


class Material(NamedObject):
    """Named global material composed from independent physical properties."""

    def __init__(
        self,
        name: str,
        *,
        elasticity: Elasticity | None = None,
        density: float | None = None,
        thermal_expansion: float | None = None,
    ) -> None:
        super().__init__(name)
        self.elasticity        = elasticity
        self.density           = None if density is None else float(density)
        self.thermal_expansion = None if thermal_expansion is None else float(thermal_expansion)

    def export(self) -> str:
        """Export the complete material definition."""

        lines = [keyword("MATERIAL", NAME=self.name)]

        if self.elasticity is not None:
            lines.append(self.elasticity.export())

        if self.density is not None:
            lines.extend([keyword("DENSITY"), csv((self.density,))])

        if self.thermal_expansion is not None:
            lines.extend([keyword("THERMALEXPANSION"), csv((self.thermal_expansion,))])

        return block(lines)


class MaterialRepository(NamedRepository[Material]):
    """Repository of named global materials."""

    def export(self) -> str:
        return "\n\n".join(material.export() for material in self)


class Profile(NamedObject):
    """General beam profile using FEMaster's scalar profile properties."""

    def __init__(
        self,
        name: str,
        area: float,
        iy: float,
        iz: float,
        j: float,
        iyz: float = 0.0,
        ey: float = 0.0,
        ez: float = 0.0,
        refy: float = 0.0,
        refz: float = 0.0,
    ) -> None:
        super().__init__(name)
        self.values = tuple(float(value) for value in (area, iy, iz, j, iyz, ey, ez, refy, refz))

    def export(self) -> str:
        return block([keyword("PROFILE", NAME=self.name), csv(self.values)])


class ProfileRepository(NamedRepository[Profile]):
    """Repository of named beam profiles."""

    def export(self) -> str:
        return "\n\n".join(profile.export() for profile in self)
