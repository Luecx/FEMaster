"""Small cantilever-style truss example using only object-valued relationships.

Every relation between FEMaster concepts is represented by the corresponding
Python object.  IDs and names are used only to create intrinsic identities and by
repositories/export; they are never passed as substitutes for nodes, regions,
materials, collectors or other referenced objects.
"""

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
    T3,
    TrussSection,
)


project = Project("cantilever")
part = project.parts.default()

# -------------------------------------------------------------------------
# Part-local topology
# -------------------------------------------------------------------------

root_node = part.nodes.add(Node(1, 0.0, 0.0, 0.0))
tip_node = part.nodes.add(Node(2, 1000.0, 0.0, 0.0))
bar_element = part.elements.add(T3(1, (root_node, tip_node)))

root = part.regions.add(NodeRegion("ROOT", (root_node,)))
tip = part.regions.add(NodeRegion("TIP", (tip_node,)))
bar = part.regions.add(ElementRegion("BAR", (bar_element,)))

# -------------------------------------------------------------------------
# Material and section
# -------------------------------------------------------------------------

steel = project.materials.add(
    Material(
        "STEEL",
        elasticity=IsotropicElasticity(210000.0, 0.3),
        density=7.85e-9,
    )
)

part.sections.add(
    TrussSection(
        "BAR_SECTION",
        element_region=bar,
        material=steel,
        area=100.0,
    )
)

# -------------------------------------------------------------------------
# Boundary conditions, load and analysis step
# -------------------------------------------------------------------------

supports = project.support_collectors.add(SupportCollector("SUPPORTS"))
supports.add(Support(root, (0.0, 0.0, 0.0)))

loads = project.load_collectors.add(LoadCollector("LOADS"))
loads.add(
    NodalForce(
        tip,
        (1000.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    )
)

project.steps.add(
    StaticStep(
        "STATIC",
        loads=(loads,),
        supports=(supports,),
    )
)

project.write("cantilever.inp")
