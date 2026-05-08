"""Minimal solid beam example."""

from __future__ import annotations

from femaster_api import FEMaster, Model
from femaster_api.model import (
    Element,
    ElementSet,
    IsotropicElasticity,
    LinearStaticLoadcase,
    LoadCollector,
    Material,
    NodalForce,
    Node,
    NodeSet,
    SolidSection,
    Support,
    SupportCollector,
)
from femaster_api.model.elements import C3D8
from femaster_api.tools.tovtk import write_vtk


if __name__ == "__main__":
    model = Model("solid_beam")

    length = 6.0
    width = 1.0
    height = 1.0

    nx = 30
    ny = 5
    nz = 5

    nodes: dict[tuple[int, int, int], Node] = {}
    for ix in range(nx + 1):
        for iy in range(ny + 1):
            for iz in range(nz + 1):
                x = length * ix / nx
                y = width * iy / ny
                z = height * iz / nz
                nodes[(ix, iy, iz)] = model.nodes.add(Node(x, y, z))

    elements = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                elements.append(
                    model.elements.add(
                        Element(
                            C3D8,
                            [
                                nodes[(ix, iy, iz)],
                                nodes[(ix + 1, iy, iz)],
                                nodes[(ix + 1, iy + 1, iz)],
                                nodes[(ix, iy + 1, iz)],
                                nodes[(ix, iy, iz + 1)],
                                nodes[(ix + 1, iy, iz + 1)],
                                nodes[(ix + 1, iy + 1, iz + 1)],
                                nodes[(ix, iy + 1, iz + 1)],
                            ],
                        )
                    )
                )

    fixed_nodes = [node for (ix, _, _), node in nodes.items() if ix == 0]
    loaded_nodes = [node for (ix, _, iz), node in nodes.items() if ix == nx and iz == nz]

    fixed = model.node_sets.add(NodeSet("FIXED", fixed_nodes))
    loaded = model.node_sets.add(NodeSet("LOADED_EDGE", loaded_nodes))
    solids = model.element_sets.add(ElementSet("BEAM_SOLIDS", elements))

    steel = model.materials.add(
        Material(
            "STEEL",
            elasticity=IsotropicElasticity(
                youngs_modulus=210000.0,
                poisson_ratio=0.3,
            ),
            density=7.85e-9,
        )
    )

    model.sections.add(
        SolidSection(
            "BEAM_SECTION",
            material=steel,
            element_set=solids,
        )
    )

    supports = model.support_collectors.add(
        SupportCollector(
            "SUPPORTS",
            [
                Support(fixed, ux=0.0, uy=0.0, uz=0.0),
            ],
        )
    )

    loads = model.load_collectors.add(
        LoadCollector(
            "LOADS",
            [
                NodalForce(loaded, fz=-500.0),
            ],
        )
    )

    model.loadcases.add(
        LinearStaticLoadcase(
            "LC_BENDING",
            supports=[supports],
            loads=[loads],
        )
    )

    result = FEMaster(model, temporary_files=True, redirect_stdout=True).run()

    print("Result summary:")
    for loadcase_name in sorted(result.loadcases):
        loadcase = result.loadcases[loadcase_name]
        print(f"  {loadcase_name}: {len(loadcase.frames)} frame(s)")
        for frame_index, frame in enumerate(loadcase.frames):
            print(f"    frame {frame_index}: {len(frame.fields)} field(s)")
            for field_name in sorted(frame.fields):
                field = frame.fields[field_name]
                preview = field.row(0) if field.data else ()
                print(f"      {field_name}: rows={field.rows}, first={preview}")

    write_vtk(model, "out.vtk", result=result)
    print("Wrote VTK: out.vtk")
