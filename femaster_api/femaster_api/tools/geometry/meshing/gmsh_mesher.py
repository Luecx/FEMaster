"""Gmsh-backed geometry mesher."""

from __future__ import annotations

from typing import Iterable

from femaster_api.model.model import Model
from femaster_api.tools.geometry.base import Curve, GeometryEntity
from femaster_api.tools.geometry.meshing.gmsh_model_builder import GmshModelBuilder
from femaster_api.tools.geometry.meshing.gmsh_to_model import append_gmsh_mesh
from femaster_api.tools.geometry.surfaces import Face
from femaster_api.tools.geometry.tagging import apply_geometry_tags
from femaster_api.tools.geometry.volumes import Volume


class GmshMesher:
    def __init__(
        self,
        *,
        mesh_size: float = 1.0,
        model_name: str = "geometry",
        recombine_transfinite: bool = True,
        fltk: bool = False,
    ) -> None:
        self.mesh_size = float(mesh_size)
        self.model_name = model_name
        self.recombine_transfinite = recombine_transfinite
        self.fltk = bool(fltk)

    def mesh(self, entities: Iterable[GeometryEntity], *, model: Model | None = None, tolerance: float = 1.0e-6) -> Model:
        gmsh = _import_gmsh()
        items = tuple(entities)
        output = Model(self.model_name) if model is None else model
        gmsh.initialize()
        try:
            gmsh.model.add(self.model_name)
            builder = GmshModelBuilder(gmsh, self.mesh_size)
            max_dim = 0
            for entity in items:
                if isinstance(entity, Face):
                    builder.add_face(entity, recombine_transfinite=self.recombine_transfinite)
                    max_dim = max(max_dim, 2)
                elif isinstance(entity, Volume):
                    builder.add_volume(entity, recombine_transfinite=self.recombine_transfinite)
                    max_dim = max(max_dim, 3)
                elif isinstance(entity, Curve):
                    builder.add_curve(entity)
                    max_dim = max(max_dim, 1)
                else:
                    raise TypeError(f"unsupported geometry entity: {type(entity).__name__}")
            if max_dim == 0:
                return output
            gmsh.model.geo.synchronize()
            gmsh.model.occ.synchronize()
            gmsh.model.mesh.generate(max_dim)
            if self.fltk:
                gmsh.fltk.run()
            append_gmsh_mesh(output, gmsh, max_dim)
        finally:
            gmsh.finalize()
        apply_geometry_tags(output, items, tolerance=tolerance)
        return output


def mesh_with_gmsh(
    entities: Iterable[GeometryEntity],
    *,
    mesh_size: float = 1.0,
    model_name: str = "geometry",
    recombine_transfinite: bool = True,
    fltk: bool = False,
    tolerance: float = 1.0e-6,
) -> Model:
    return GmshMesher(
        mesh_size=mesh_size,
        model_name=model_name,
        recombine_transfinite=recombine_transfinite,
        fltk=fltk,
    ).mesh(
        entities,
        tolerance=tolerance,
    )


def _import_gmsh():
    try:
        import gmsh
    except ImportError as exc:
        raise RuntimeError("gmsh is required for geometry meshing") from exc
    return gmsh
