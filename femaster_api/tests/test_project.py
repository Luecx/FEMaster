"""Behavioral checks for project ownership, export and native INP reading."""

from femaster_api import (
    ElementRegion,
    ElementSurface,
    IsotropicElasticity,
    LoadCollector,
    Material,
    Node,
    NodeRegion,
    NodalForce,
    Project,
    StaticStep,
    Support,
    SupportCollector,
    SurfaceRegion,
    T3D2,
    TrussSection,
)


def make_project() -> Project:
    project = Project("TEST")
    part = project.parts.default()

    part.nodes.add(Node(10, 0.0, 0.0, 0.0))
    part.nodes.add(Node(20, 1000.0, 0.0, 0.0))
    part.elements.add(T3D2(50, (10, 20)))

    part.regions.add(NodeRegion("ROOT", (10,)))
    part.regions.add(NodeRegion("TIP", (20,)))
    part.regions.add(ElementRegion("BAR", (50,)))
    part.regions.add(SurfaceRegion("SURFACE_SET", ("OUTER",)))
    part.surfaces.add(ElementSurface("OUTER").add(50, 1))

    project.materials.add(
        Material(
            "STEEL",
            elasticity=IsotropicElasticity(210000.0, 0.3),
        )
    )
    part.sections.add(
        TrussSection("BAR_SECTION", "BAR", "STEEL", 100.0)
    )

    supports = project.support_collectors.add(SupportCollector("BC"))
    supports.add(Support("ROOT", (0.0, 0.0, 0.0)))

    loads = project.load_collectors.add(LoadCollector("LOAD"))
    loads.add(
        NodalForce(
            "TIP",
            (1000.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        )
    )

    project.steps.add(
        StaticStep(
            "STATIC",
            loads=("LOAD",),
            supports=("BC",),
        )
    )
    return project


def test_default_part_is_permanent_repository_index_zero():
    project = Project()

    assert project.parts.default() is project.parts[0]

    try:
        del project.parts[0]
    except ValueError:
        pass
    else:
        raise AssertionError("default part must not be removable")


def test_project_export_and_read_inp_text_round_trip_core_semantics():
    source = make_project()
    text = source.export()

    assert "*SURFACE, NAME=OUTER, TYPE=ELEMENT" in text
    assert "*SFSET, SFSET=SURFACE_SET" in text
    assert "*LOADCASE, TYPE=LINEARSTATIC, NAME=STATIC" in text

    project = Project.read_inp_text(text)
    part = project.parts.default()

    assert project.name == "TEST"
    assert part.nodes[10].x == 0.0
    assert part.nodes[20].x == 1000.0
    assert part.elements[50].nodes == (10, 20)
    assert part.surfaces["OUTER"].entries == [(50, 1)]
    assert part.regions.surfaces["SURFACE_SET"].members == ["OUTER"]
    assert project.materials["STEEL"].elasticity.youngs_modulus == 210000.0
    assert project.steps["STATIC"].loads == ("LOAD",)
