# Bolt pretension

The complete keyword reference is in the User Manual chapter **Bolt Pretension**
(`documentation/pages/11-pretension.tex`).

Pretension supports linear-static analysis of planar, matching solid-element
faces. A merged mesh is separated along existing element boundaries; elements
are not cut through their interiors. Keep supports away from the cut.

## Native input

Define a cut that completely separates two solid regions:

```text
*SURFACE, SFSET=CUT_FACE
2, 2, 1
*PRETENSIONSECTION, NAME=PT1, SURFACE_A=CUT_FACE
```

The surface row is `surface-id, element-id, local-face-id`. For two already
separated, coincident faces, add `SURFACE_B=CUT_B` and use `SURFACE_A=CUT_A`.
Native sections generate an isolated scalar control node automatically.

Inside the first linear-static loadcase, apply the total section force:

```text
*PRETENSION, SECTION=PT1, ACTION=LOAD, CONTROL=FORCE, VALUE=10000
*CONSTRAINTMETHOD, TYPE=ELIMINATION
```

Alternatively prescribe the axial adjustment with
`CONTROL=DISPLACEMENT, VALUE=0.001`. In the next loadcase use:

```text
*PRETENSION, SECTION=PT1, ACTION=LOCK
```

LOCK holds the solved adjustment, not the section force. Later external loads
can change the transmitted force while the adjustment remains fixed. Native
supports and named external load collectors must be selected in the applicable
loadcases. A CLOAD on a node set is a force per node; PRETENSION VALUE is the
entire section force. Native SUPPORT rows contain displacement components,
not Abaqus DOF ranges.

The axis follows the first SURFACE_A face normal. The relative displacement is
`q = (u_B - u_A) dot axis`; positive force is work-conjugate to positive q.
Transverse relative displacements remain tied. Use NULLSPACE or ELIMINATION;
LAGRANGE is rejected when a pretension section is present.

## Abaqus-compatible input

Declare an isolated control node explicitly and reference the cut surface:

```text
*PRE-TENSION SECTION, SURFACE=CUT, NODE=900001, FOLLOWER=NO
```

In the first linear-static STEP, `*CLOAD` with `900001, 1, 10000` applies the
force; `*BOUNDARY` with `900001, 1, 1, 0.001` prescribes the adjustment.
In a later step, `*BOUNDARY, FIXED` with `900001, 1, 1` holds the preceding
successful adjustment. A carried CLOAD remains present with FIXED. Native
PRETENSION actions and Abaqus control conditions cannot be mixed on one section.
The manual documents history replacement, scopes, normal overrides and limits.

## Results and examples

PTGAP measures the axial gap; PTFORC contains interface forces and the signed
section resultant in control-node component 1; TF stores the total resultant
at the control node. Results are available in RES, FRD and FEMR.

- [Native force / LOCK / external-load example](../examples/09_constraints/095_bolt_pretension.inp)
- [Abaqus-compatible companion](../examples/09_constraints/096_bolt_pretension_abaqus.inp)
- [Validation scope and commands](pretension-validation.md)
- [Reproducible preparation and solver benchmarks](pretension-benchmark.md)
