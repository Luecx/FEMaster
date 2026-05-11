from femaster_api import *

model = Model("cantilever")

# ------------------------------------------------------------------
# Nodes
# ------------------------------------------------------------------

n1 = model.nodes.add(Node(0.0, 0.0, 0.0))
n2 = model.nodes.add(Node(1.0, 0.0, 0.0))

# ------------------------------------------------------------------
# Element
# ------------------------------------------------------------------

beam = model.elements.add(Element(B33, (n1, n2)))

# ------------------------------------------------------------------
# Sets
# ------------------------------------------------------------------

beam_set  = model.sets.add(ElementSet("BEAM", (beam,)))
fixed_set = model.sets.add(NodeSet("FIXED", (n1,)))
tip_set   = model.sets.add(NodeSet("TIP", (n2,)))

# ------------------------------------------------------------------
# Material
# ------------------------------------------------------------------

steel = Material("STEEL")
steel = steel.set_isotropic_elasticity(
    E  = 210e9,
    nu = 0.3,
)
steel = steel.set_density(7850)

steel = model.materials.add(steel)

# ------------------------------------------------------------------
# Beam profile
# ------------------------------------------------------------------

profile = Profile(
    name = "BOX",
    area = 2.0e-4,
    iy   = 1.0e-8,
    iz   = 2.0e-8,
    j    = 5.0e-9,
)

profile = model.profiles.add(profile)

# ------------------------------------------------------------------
# Section
# ------------------------------------------------------------------

section = BeamSection(
    name        = "BEAM_SECTION",
    material    = steel,
    element_set = beam_set,
    profile     = profile,
    orientation = (0.0, 1.0, 0.0),
)

model.sections.add(section)

# ------------------------------------------------------------------
# Boundary conditions
# ------------------------------------------------------------------

bc = SupportCollector("BCS").add(
    Support(
        target = fixed_set,
        values = (0, 0, 0, 0, 0, 0),
    )
)

model.support_collectors.add(bc)

# ------------------------------------------------------------------
# Loads
# ------------------------------------------------------------------

loads = LoadCollector("LOADS").add(
    NodalForce(
        target = tip_set,
        values = (0.0, 0.0, -1000.0, 0.0, 0.0, 0.0),
    )
)

model.load_collectors.add(loads)

# ------------------------------------------------------------------
# Analysis step
# ------------------------------------------------------------------

step = StaticStep(
    name     = "STATIC",
    loads    = (loads,),
    supports = (bc,),
)

model.steps.add(step)

# ------------------------------------------------------------------
# Validate
# ------------------------------------------------------------------

model.validate(raise_on_error=True)

# ------------------------------------------------------------------
# Export
# ------------------------------------------------------------------

results = FEMaster(model=model).run("cantilever.inp")

displacement = results.field("DISPLACEMENT")

print(displacement.row(0))
