"""FEMaster serialization for non-element model features."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.features import FeatureRepository


def write_features(features: FeatureRepository) -> str:
    blocks: list[str] = []
    for point_mass in features.point_masses():
        blocks.append(
            block(
                [
                    keyword("POINTMASS", NSET=point_mass.node_set.name),
                    csv((point_mass.mass, *point_mass.inertia, *point_mass.spring, *point_mass.rotational_spring)),
                ]
            )
        )
    return "\n\n".join(block for block in blocks if block)
