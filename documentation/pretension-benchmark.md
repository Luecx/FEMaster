# Pretension benchmark protocol

The benchmarks separate mesh preparation, node matching and equilibrium solving.
Timings depend on hardware, build flags, solver backend and background load.
They are descriptive measurements, not correctness thresholds.

## Solver-free preparation

Configure with `FEMASTER_BUILD_BENCHMARKS=ON`, then build:

```sh
cmake --build build/pretension-bench --config Release --target FEMaster_pretension_benchmark FEMaster_pretension_matching
```

The executables are in `bin` for single-config generators or `bin/Release` for
Visual Studio. The preparation benchmark builds fresh structured C3D8 bars,
three elements long, with n by n elements across and one planar merged cut.
It times only `prepare_pretension_sections()`: compilation, material assignment,
stiffness assembly, equilibrium solving and result writing are excluded.
Preparation diagnostic logging is included.

Sizes are n=4,8,16,32,64. Repeat 0 is warmup; repeats 1-5 are measured samples.
Report the median and range for each size and preserve the raw output with the
exact source revision and build settings. The generated control node is included
in final node counts.

The matching executable measures an isolated copy of the deterministic matching
kernel. It excludes topology work and is a diagnostic lower-level measurement,
not an end-to-end preparation benchmark.

## Equilibrium and result comparison

Generate the corresponding analytical bar models:

```sh
python benchmarks/pretension_solver.py build/pretension-models
python benchmarks/pretension_compare.py BASELINE_EXE CANDIDATE_EXE build/pretension-models build/pretension-comparison.csv
```

Replace BASELINE_EXE and CANDIDATE_EXE with two actual executable paths. Both
must support the native pretension input and use the same precision, backend
and threading settings. The comparison script records the CLI's separate
`solving constraint system` timer, checks analytical gap and force, and compares
exported RES data. It runs the baseline and candidate sequentially, with one
warmup and five samples per size. Use an otherwise idle machine.

Record compiler/version, optimization flags, Eigen/MKL versions, constraint
method, thread count, source revisions and full command lines with each result.
Do not interpret a change of sparse solver or thread count as a preparation
speedup. A total-process timing also includes parsing, constraints, recovery and
output. Historical development measurements are not required to reproduce the
current protocol.
