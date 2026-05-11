"""Create a small geometry mesh with optional Gmsh FLTK display."""

from __future__ import annotations

import argparse

from femaster_api.model.sets import EntityType
from femaster_api.tools.geometry import Face, Line, Volume, mesh_with_gmsh
from femaster_api.tools import tovtk


def main(*, fltk: bool = False) -> None:
    volume = _cube(tags=("SOLID",))
    model = mesh_with_gmsh(
        (volume,),
        mesh_size=0.75,
        model_name="geometry_meshing_example",
        recombine_transfinite=False,
        fltk=fltk,
    )

    solid_nodes = model.sets.get(EntityType.NODE, "SOLID")
    solid_elements = model.sets.get(EntityType.ELEMENT, "SOLID")
    print(f"nodes: {len(model.nodes)}")
    print(f"elements: {len(model.elements)}")
    print(f"SOLID node set: {len(solid_nodes.members)}")
    print(f"SOLID element set: {len(solid_elements.members)}")

    tovtk(model, "geometry_meshing_example.vtk")



def _cube(*, tags: tuple[str, ...] = ()) -> Volume:
    p000 = (0.0, 0.0, 0.0)
    p100 = (1.0, 0.0, 0.0)
    p110 = (1.0, 1.0, 0.0)
    p010 = (0.0, 1.0, 0.0)
    p001 = (0.0, 0.0, 1.0)
    p101 = (1.0, 0.0, 1.0)
    p111 = (1.0, 1.0, 1.0)
    p011 = (0.0, 1.0, 1.0)
    return Volume(
        (
            _quad(p000, p100, p110, p010),
            _quad(p001, p011, p111, p101),
            _quad(p000, p001, p101, p100),
            _quad(p010, p110, p111, p011),
            _quad(p000, p010, p011, p001),
            _quad(p100, p101, p111, p110),
        ),
        tags=tags,
    )


def _quad(a, b, c, d) -> Face:
    return Face(
        (
            Line(a, b, seed=12),
            Line(b, c, seed=12),
            Line(c, d, seed=12),
            Line(d, a, seed=12),
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fltk", action="store_true")
    args = parser.parse_args()
    main(fltk=args.fltk)
