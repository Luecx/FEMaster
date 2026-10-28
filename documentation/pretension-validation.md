# Pretension validation scope

## Verified review state (27 October 2026)

The focused target `FEMaster_pretension_test` built successfully with MSVC
19.44.35228, Release, double precision, CPU DIRECT, OpenMP and MKL disabled.
All **27 tests passed** after integration into the common condition history.
The core and focused target were rebuilt from clean objects after the
ConditionState layout change. The normal FEMaster executable built successfully,
and the CLI face-pair loadstep regression passed as well. The tested source copy
matched the workspace's `src` and `tests` files (ignoring Windows line endings).

The full `FEMaster_test` target remains blocked at compilation. Comparing the
merge commit `00b31a55` with its unmodified upstream parent `06b15a0e` reproduced
all 34 distinct rework test compiler diagnostics in upstream: outdated element
and result-recovery calls, plus `rand_r` under MSVC. Upstream additionally failed
in `element_solid.ipp:285` with C2244; that signature correction is already in
the rework branch. The broad suite was not run or rebuilt again after the
condition-history refactor; the earlier compiler failures remain known blockers.
These focused results do not prove absence of other integration defects.

## Analytical and negative regressions

The linear static regression uses two coincident faces with separate node IDs.
For the straight bar with E=210000, nu=0, A=4 and effective length 3, a total
force of 1000 gives an expected axial gap of 0.0035714285714285714.

[The focused test](../tests/test_pretension.cpp) checks FORCE, HOLD, external
loading, common axial displacement, no added nodes, transmitted section force,
and rejection of shared topology, missing matches, ambiguous matches and
HOLD before a successfully completed step.
Negative post-check tests also verify that a rejected first step cannot enable
HOLD, while a rejected later step leaves the previously accepted gap unchanged.
Condition-history regressions check frozen initial prescriptions, replacement
without accumulation, inherited FORCE and HOLD, and a rejected target that must
not become the next step's inherited prescription. FORCE/HOLD/FORCE cycles are
checked through both LOADCASE and STEP/STATIC input syntax.
[The CLI regression](../tests/run_pretension_face_pair.cmake) checks all three
steps against [the input deck](../tests/data/pretension_face_pair_loadsteps.inp).

The focused suite also checks rotated geometry and loads with NULLSPACE and
ELIMINATION. Every pair must retain the common axial gap and zero tangential
relative displacement in FORCE and HOLD; PTFORC must be equal and opposite,
and its section sum must agree with TF and the analytical force.
For a total axial load of -4 at length coordinate 2, the held cut at coordinate 1
transmits 1000 + 4/3; doubling that load gives 1000 + 8/3.
The [two-section deck](../tests/data/pretension_two_sections.inp) uses independent
orthogonal bars with pretension forces 1000 and 2000. Loading the first bar must
not change the second bar's force or gap. Both constraint methods are covered.

The implementation is limited to planar, coincident solid faces and linear
static steps. Additional active constraints on either cut face are rejected;
external loads remain supported. Nonlinear, perturbation and LAGRANGE pretension
are unsupported.
No external solver comparison has been completed for this revision.
