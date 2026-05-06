"""FEMaster serialization for surfaces."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword, target_token
from femaster_api.model.surfaces import SurfaceRepository


def write_surfaces(surfaces: SurfaceRepository) -> str:
    blocks: list[str] = []
    for surface_set in surfaces.sets():
        lines = [keyword("SURFACE", NAME=surface_set)]
        for surface in surfaces.by_set(surface_set):
            if surface.id is None:
                lines.append(csv((target_token(surface.target), surface.side)))
            else:
                lines.append(csv((surface.id, target_token(surface.target), surface.side)))
        blocks.append(block(lines))
    return "\n\n".join(block for block in blocks if block)
