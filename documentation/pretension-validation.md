# Pretension validation scope

The feature suite is independently runnable when `FEMASTER_BUILD_TESTS=ON`:

```sh
cmake --build build/pretension-tests --config Release --target FEMaster FEMaster_pretension_test
ctest --test-dir build/pretension-tests -C Release -R Pretension --output-on-failure
```

Configure that build directory first using the platform's normal compiler and
dependency settings. For single-config generators, configuration is selected
at configure time. No external solver is required for these tests.

## Automated evidence

- Analytical straight bar with E=210000, nu=0, A=4, L=3: a 1000 force gives
  `q = F L / (E A) = 0.0035714285714285714`.
- Force and adjustment control, successful-solve LOCK and subsequent external
  load with NULLSPACE and ELIMINATION.
- Native/Abaqus-compatible input equivalence; carried force, FIXED, OP=MOD/NEW,
  amplitudes and scalar node/singleton alias replacement.
- Mixed and redundant constraints; rejection of contradictory prescriptions,
  unsupported analysis/constraint methods, invalid controls and invalid cuts.
- Two independent sections, rotated cuts, generated/non-dense result IDs,
  Assembly/Instance references and owned boundary rewiring.
- RES/FRD/FEMR field selection, control-node accounting and binary/RES values.

These tests establish the stated formulation and supported input paths. They
do not establish complete Abaqus or CalculiX compatibility.

## Upstream regression comparison

The Windows/MSVC Release comparison used upstream
`22cb7c8b2ec56fee9cea4622ecfaa02accb17a75`, double precision, identical dependency
revisions and identical compiler flags. Both runnable general test subsets
contained 770 tests: 752 passed and the same 18 failed. The complete general
test target could not be built on either revision: eight test files contain
obsolete operator calls or Windows-incompatible `rand_r` use. Those files were
not silently repaired or counted as passing.

Thirteen common CLI examples were executed. Twelve produced numerically
identical serialized RES fields; four of those already reported failed
post-checks on upstream. The thirteenth was rejected by both programs because
nonlinear follower pressure is unsupported. This is bounded regression
evidence, not a claim that the entire upstream suite is green.

A separate MSVC Mat6 declaration/definition compatibility adjustment was needed
locally to build the reference. It is outside the proposed pretension patch.
The C3D15 S2 and C3D20 S6 surface-connectivity corrections and their general
topology test are likewise separate general-element fixes. Validation models
using those faces require the corresponding correction.

## Clean contribution build

The contribution branch also builds on Windows with IntelLLVM 2026.1.1,
Ninja, Release, double precision, MKL disabled and OpenMP disabled. The CLI
and all 39 pretension tests pass without the separate Mat6 or surface fixes.
The benchmark executables and the generated straight-bar solver case also run
successfully. Activate the Visual Studio and oneAPI environment before building
and running tests so that the Intel runtime libraries are available.

## Independent solver comparison: pending

A manually exercised screw model has plausible force/preload, locked-adjustment
and service-load steps. Its displayed stress and physical displacement maxima
agree closely with CalculiX in the inspected force case. A reproducible numeric
comparison of physical displacements, section resultants and support reactions
for all three steps, plus a separate prescribed-adjustment case, remains to be
supplied by the model owner. No completed external validation is claimed here.

Runtime tuning and automatic deformation scaling in the visualizer are outside
this contribution. Reproducible benchmark commands are documented separately
in [pretension-benchmark.md](pretension-benchmark.md).
