"""FEMaster serialization for surfaces."""

from __future__ import annotations

from femaster_api.export.context import ExportContext
from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.surfaces import SurfaceRepository


def write_surfaces(surfaces: SurfaceRepository, context: ExportContext) -> str:
    blocks: list[str] = []
    for surface_set in surfaces.sets():
        lines = [keyword("SURFACE", NAME=surface_set)]
        for surface in surfaces.by_set(surface_set):
            lines.append(csv((context.target_token(surface.target), surface.side)))
        blocks.append(block(lines))
    return "\n\n".join(block for block in blocks if block)
