"""Behavioral checks for object-valued model relationships and INP round trips."""

import pytest

from femaster_api import (
    Amplitude,
    ElementRegion,
    ElementSurface,
    IsotropicElasticity,
    LoadCollector,
    Material,
    Node,
    NodeRegion,
    NodalForce,
    PressureLoad,
    Project,
    RectangularCoordinateSystem,
    StaticStep,
    Support,
    SupportCollector,
    SurfaceRegion,
    T3,
    TrussSection,
)


def make_project() -> Project:
    project = Project("TEST")
    part = project.parts.default()

    node_root = part.nodes.add(Node(10, 0.0, 0.0, 0.0))
    node_tip = part.nodes.add(Node(20, 1000.0, 0.0, 0.0))
    element = part.elements.add(T3(50, (node_root, node_tip)))

    root = part.regions.add(NodeRegion("ROOT", (node_root,)))
    tip = part.regions.add(NodeRegion("TIP", (node_tip,)))
    bar = part.regions.add(ElementRegion("BAR", (element,)))
    outer = part.surfaces.add(ElementSurface("OUTER").add(element, 1))
    outer_group = part.regions.add(SurfaceRegion("SURFACE_SET", (outer,)))

    material = project.materials.add(
        Material(
            "STEEL",
            elasticity=IsotropicElasticity(210000.0, 0.3),
        )
    )
    part.sections.add(
        TrussSection("BAR_SECTION", bar, material, 100.0)
    )

    orientation = project.coordinate_systems.add(
        RectangularCoordinateSystem(
            "GLOBAL",
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
        )
    )
    amplitude = project.amplitudes.add(
        Amplitude("RAMP").add(0.0, 0.0).add(1.0, 1.0)
    )

    supports = project.support_collectors.add(SupportCollector("BC"))
    supports.add(Support(root, (0.0, 0.0, 0.0), orientation=orientation))

    loads = project.load_collectors.add(LoadCollector("LOAD"))
    loads.add(
        NodalForce(
            tip,
            (1000.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            orientation=orientation,
            amplitude=amplitude,
        )
    )
    loads.add(PressureLoad(outer_group, 2.0, amplitude=amplitude))

    project.steps.add(
        StaticStep(
            "STATIC",
            loads=(loads,),
            supports=(supports,),
        )
    )
    return project


def test_default_part_is_permanent_repository_index_zero():
    project = Project()

    assert project.parts.default() is project.parts[0]

    with pytest.raises(ValueError):
        del project.parts[0]


def test_project_export_and_read_inp_text_round_trip_object_relationships():
    source = make_project()
    text = source.export()

    assert "*ELEMENT, TYPE=T3" in text
    assert "*SURFACE, NAME=OUTER, TYPE=ELEMENT" in text
    assert "*SFSET, SFSET=SURFACE_SET" in text
    assert "AMPLITUDE=RAMP" in text
    assert "ORIENTATION=GLOBAL" in text
    assert "*LOADCASE, TYPE=LINEARSTATIC, NAME=STATIC" in text

    project = Project.read_inp_text(text)
    part = project.parts.default()

    assert project.name == "TEST"
    assert part.elements[50].nodes == (part.nodes[10], part.nodes[20])
    assert part.regions.nodes["ROOT"].members == [part.nodes[10]]
    assert part.regions.elements["BAR"].members == [part.elements[50]]
    assert part.surfaces["OUTER"].entries == [(part.elements[50], 1)]
    assert part.regions.surfaces["SURFACE_SET"].members == [part.surfaces["OUTER"]]

    section = part.sections["BAR_SECTION"]
    assert section.element_region is part.regions.elements["BAR"]
    assert section.material is project.materials["STEEL"]

    support = project.support_collectors["BC"].supports[0]
    assert support.target is part.regions.nodes["ROOT"]
    assert support.orientation is project.coordinate_systems["GLOBAL"]

    force = project.load_collectors["LOAD"].loads[0]
    assert force.target is part.regions.nodes["TIP"]
    assert force.orientation is project.coordinate_systems["GLOBAL"]
    assert force.amplitude is project.amplitudes["RAMP"]

    pressure = project.load_collectors["LOAD"].loads[1]
    assert pressure.target is part.regions.surfaces["SURFACE_SET"]
    assert pressure.amplitude is project.amplitudes["RAMP"]

    step = project.steps["STATIC"]
    assert step.loads == (project.load_collectors["LOAD"],)
    assert step.supports == (project.support_collectors["BC"],)


def test_cross_object_relationships_reject_string_surrogates():
    project = make_project()
    part = project.parts.default()
    node = part.nodes[10]
    amplitude = project.amplitudes["RAMP"]

    with pytest.raises(TypeError):
        NodalForce("TIP")

    with pytest.raises(TypeError):
        NodalForce(node, amplitude="RAMP")

    with pytest.raises(TypeError):
        T3(99, (10, 20))

    with pytest.raises(TypeError):
        NodeRegion("BAD", ("10",))

    with pytest.raises(TypeError):
        StaticStep("BAD", loads=("LOAD",))

    valid = NodalForce(node, amplitude=amplitude)
    assert valid.target is node
    assert valid.amplitude is amplitude
