"""Behavioral checks for object-valued model relationships and INP round trips."""

import pytest

from femaster_api import (
    Amplitude,
    ElementRegion,
    ElementSurface,
    Field,
    FieldDomain,
    Instance,
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

    density = project.fields.add(
        Field("RHO", FieldDomain.ELEMENT, ("RHO",))
    )
    density.set(element, (0.8,))

    directors = project.fields.add(
        Field("DIRECTORS", FieldDomain.ELEMENT_NODAL, ("X", "Y", "Z"))
    )
    directors.set((element, 0), (0.0, 0.0, 1.0))
    directors.set((element, 1), (0.0, 0.0, 1.0))

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
    assert "*FIELD, NAME=RHO, TYPE=ELEMENT" in text
    assert "50, 0.8" in text

    project = Project.read_inp_text(text)
    part = project.parts.default()

    assert project.name == "TEST"
    assert part.elements[50].nodes == (part.nodes[10], part.nodes[20])
    assert part.regions.nodes["ROOT"].members == [part.nodes[10]]
    assert part.regions.elements["BAR"].members == [part.elements[50]]
    assert part.surfaces["OUTER"].entries == [(part.elements[50], 1)]
    assert part.regions.surfaces["SURFACE_SET"].members == [part.surfaces["OUTER"]]

    # Native section keywords identify assignments by ELSET and do not persist
    # the Python repository name, so validate the reconstructed section by order
    # and, importantly, by object identity of its semantic relationships.
    section = part.sections[0]
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

    density = project.fields["RHO"]
    assert list(density.values) == [part.elements[50]]
    assert density[part.elements[50]] == (0.8,)

    directors = project.fields["DIRECTORS"]
    assert set(directors.values) == {
        (part.elements[50], 0),
        (part.elements[50], 1),
    }


def test_cross_object_relationships_reject_string_and_id_surrogates():
    project = make_project()
    part = project.parts.default()
    node = part.nodes[10]
    element = part.elements[50]
    material = project.materials["STEEL"]
    amplitude = project.amplitudes["RAMP"]

    with pytest.raises(TypeError):
        NodalForce("TIP")

    with pytest.raises(TypeError):
        NodalForce(node, amplitude="RAMP")

    with pytest.raises(TypeError):
        T3(99, (10, 20))

    with pytest.raises(TypeError):
        NodeRegion("BAD", (10,))

    with pytest.raises(TypeError):
        ElementRegion("BAD", (50,))

    with pytest.raises(TypeError):
        TrussSection("BAD", "BAR", material, 1.0)

    with pytest.raises(TypeError):
        Support("ROOT", (0.0,))

    with pytest.raises(TypeError):
        PressureLoad("OUTER", 1.0)

    with pytest.raises(TypeError):
        StaticStep("BAD", loads=("LOAD",))

    with pytest.raises(TypeError):
        StaticStep("BAD", supports=("BC",))

    with pytest.raises(TypeError):
        Instance("BAD", "PART")

    element_field = Field("RHO2", FieldDomain.ELEMENT, ("RHO",))
    with pytest.raises(TypeError):
        element_field.set(element.id, (1.0,))

    nodal_field = Field("TEMP", FieldDomain.NODE, ("T",))
    with pytest.raises(TypeError):
        nodal_field.set(node.id, (20.0,))

    local_field = Field("IP", FieldDomain.ELEMENT_IP, ("V",))
    with pytest.raises(TypeError):
        local_field.set((element.id, 0), (1.0,))

    valid = NodalForce(node, amplitude=amplitude)
    assert valid.target is node
    assert valid.amplitude is amplitude


def test_model_field_local_indices_remain_intrinsic_integers():
    project = make_project()
    element = project.parts.default().elements[50]

    field = Field("MP", FieldDomain.ELEMENT_MP, ("V",))
    field.set((element, 2, 3), (7.0,))

    assert field[(element, 2, 3)] == (7.0,)
    assert "50, 2, 3, 7.0" in field.export()

    with pytest.raises(ValueError):
        field.set((element, -1, 0), (1.0,))
