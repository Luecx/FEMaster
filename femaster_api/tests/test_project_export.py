from femaster_api import (
    ElementRegion,
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


def test_project_export_is_femaster_deck():
    text = make_project().export()

    assert "*MODEL, NAME=TEST" in text
    assert "*NODE" in text
    assert "*ELEMENT, TYPE=T3D2" in text
    assert "*MATERIAL, NAME=STEEL" in text
    assert "*TRUSSSECTION, ELSET=BAR, MATERIAL=STEEL" in text
    assert "*SUPPORT, SUPPORT_COLLECTOR=BC" in text
    assert "*CLOAD, LOAD_COLLECTOR=LOAD" in text
    assert "*LOADCASE, TYPE=LINEARSTATIC, NAME=STATIC" in text
    assert text.rstrip().endswith("*END")
