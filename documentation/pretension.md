# Bolt pretension

Bolt pretension is supported in linear static loadcases. Define two coincident,
planar solid faces with separate node IDs. The preprocessor must split the mesh;
FEMaster does not change topology. Every node on one face must have exactly one
node on the other face at the same reference coordinates.

```text
*SURFACE, SFSET=CUT_A
1, 1, 2
*SURFACE, SFSET=CUT_B
2, 2, 1
*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_A, SURFACE_B=CUT_B
```

The axis follows the area-weighted normal of SURFACE_A. Both tangential
components of relative displacement are fixed for each pair. In FORCE mode,
all pairs share one free axial relative displacement. The specified total
force is distributed equally across the pairs, with opposite forces on the
two faces. No control node or additional degree of freedom is created.

```text
*LOADCASE, TYPE=LINEARSTATIC, NAME=PRELOAD
*PRETENSION, SECTION=PT1, MODE=FORCE, VALUE=1000
*END
*LOADCASE, TYPE=LINEARSTATIC, NAME=HOLD
*PRETENSION, SECTION=PT1, MODE=HOLD
*END
```

HOLD requires a successfully completed preceding static step. It fixes the
axial gap reached in that step; subsequent external loads may change the
transmitted force. A section with no FORCE or HOLD action starts with zero gap.
Use NULLSPACE or ELIMINATION constraint handling. LAGRANGE, nonlinear and
perturbation steps with pretension are not supported.

Failed numerical post-checks reject the step before results and step history
are accepted. The last accepted pretension gap remains available to HOLD;
HOLD without a previously accepted gap is rejected.

FORCE and HOLD are immutable conditions in the same initial/target step-history
framework as other prescriptions. A new command replaces the target condition
for its section without changing the frozen initial condition. Omitting a new
command inherits the previously accepted mode and value. HOLD freezes the last
accepted gap when its condition is defined. A rejected step does not become the
inherited input state for the next step. Section geometry is independent of
these prescriptions; only accepted solution results are stored on the section.

PTGAP gives the axial relative displacement on both nodes of each pair.
PTFORC gives equal and opposite transmitted force vectors on the two faces.
The signed section force is the sum of PTFORC on side B projected onto the
section axis. TF repeats that total at the side B nodes for convenient
inspection. Active supports, couplings, MPCs, ties, connectors, rigid-body constraints and manual equations
on either cut face are rejected because their forces cannot be separated from the interface force.
Place these constraints away from the cut. External loads at cut nodes remain supported.

See [the complete example](../examples/09_constraints/095_bolt_pretension.inp).
