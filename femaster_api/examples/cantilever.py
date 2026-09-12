"""Small cantilever-style truss example using the rewritten FEMaster API."""

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


project = Project("cantilever")
part = project.parts.default()

# -------------------------------------------------------------------------
# Part-local topology
# -------------------------------------------------------------------------

part.nodes.add(Node(1, 0.0, 0.0, 0.0))
part.nodes.add(Node(2, 1000.0, 0.0, 0.0))
part.elements.add(T3D2(1, (1, 2)))

part.regions.add(NodeRegion("ROOT", (1,)))
part.regions.add(NodeRegion("TIP", (2,)))
part.regions.add(ElementRegion("BAR", (1,)))

# -------------------------------------------------------------------------
# Material and section
# -------------------------------------------------------------------------

project.materials.add(
    Material(
        "STEEL",
        elasticity=IsotropicElasticity(210000.0, 0.3),
        density=7.85e-9,
    )
)

part.sections.add(
    TrussSection(
        "BAR_SECTION",
        element_region="BAR",
        material="STEEL",
        area=100.0,
    )
)

# -------------------------------------------------------------------------
# Boundary conditions and load
# -------------------------------------------------------------------------

supports = project.support_collectors.add(SupportCollector("SUPPORTS"))
supports.add(Support("ROOT", (0.0, 0.0, 0.0)))

loads = project.load_collectors.add(LoadCollector("LOADS"))
loads.add(NodalForce("TIP", (1000.0, 0.0, 0.0, 0.0, 0.0, 0.0)))

project.steps.add(
    StaticStep(
        "STATIC",
        loads=("LOADS",),
        supports=("SUPPORTS",),
    )
)

project.write("cantilever.inp")
