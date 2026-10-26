# Pretension validation scope

The linear static regression uses two coincident faces with separate node IDs.
For the straight bar with E=210000, nu=0, A=4 and effective length 3, a total
force of 1000 gives an expected axial gap of 0.0035714285714285714.

[The focused test](../tests/test_pretension.cpp) checks FORCE, HOLD, external
loading, common axial displacement, no added nodes, transmitted section force,
and rejection of shared topology, missing matches, ambiguous matches and
HOLD before a successfully completed step.
[The CLI regression](../tests/run_pretension_face_pair.cmake) checks all three
steps against [the input deck](../tests/data/pretension_face_pair_loadsteps.inp).

The implementation is limited to planar, coincident solid faces and linear
static steps. Supports on the cut nodes make the reported interface force
ambiguous. Nonlinear, perturbation and LAGRANGE pretension are unsupported.
No external solver comparison has been completed for this revision.