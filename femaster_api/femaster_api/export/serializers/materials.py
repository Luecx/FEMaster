"""FEMaster serialization for materials."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.materials import (
    ABDElasticity,
    GeneralizedIsotropicElasticity,
    IsotropicElasticity,
    MaterialRepository,
    OrthotropicElasticity,
)


def write_materials(materials: MaterialRepository) -> str:
    blocks = [_write_material(material) for material in materials.all()]
    return "\n\n".join(block for block in blocks if block)


def _write_material(material) -> str:
    lines = [keyword("MATERIAL", NAME=material.name)]
    elasticity = material.elasticity
    if isinstance(elasticity, IsotropicElasticity):
        lines.append(keyword("ELASTIC", TYPE="ISOTROPIC"))
        lines.append(csv((elasticity.youngs_modulus, elasticity.poisson_ratio)))
    elif isinstance(elasticity, GeneralizedIsotropicElasticity):
        lines.append(keyword("ELASTIC", TYPE="GENISO"))
        lines.append(csv((elasticity.youngs_modulus, elasticity.poisson_ratio, elasticity.shear_modulus)))
    elif isinstance(elasticity, OrthotropicElasticity):
        lines.append(keyword("ELASTIC", TYPE="ORTHOTROPIC"))
        lines.append(
            csv(
                (
                    elasticity.e1,
                    elasticity.e2,
                    elasticity.e3,
                    elasticity.g23,
                    elasticity.g13,
                    elasticity.g12,
                    elasticity.nu23,
                    elasticity.nu13,
                    elasticity.nu12,
                )
            )
        )
    elif isinstance(elasticity, ABDElasticity):
        lines.append(keyword("ELASTIC", TYPE="ABD"))
        values = elasticity.values
        for start in range(0, len(values), 8):
            lines.append(csv(values[start : start + 8]))

    if material.density is not None:
        lines.append(keyword("DENSITY"))
        lines.append(csv((material.density,)))
    if material.thermal_expansion is not None:
        lines.append(keyword("THERMALEXPANSION"))
        lines.append(csv((material.thermal_expansion,)))
    return block(lines)
