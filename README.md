# FEMaster

FEMaster is a C++17 finite-element solver for structural mechanics. One CMake
build now covers Linux/WSL, Windows, and macOS. Intel oneMKL is supported on
Linux and Windows; CUDA and optional cuDSS are supported on Linux and Windows.

## Features

- Linear and geometrically nonlinear static analysis
- Eigenfrequency, buckling, and implicit Newmark analyses
- Beam, shell, solid, truss, and point-mass models
- Supports, couplings, connectors, ties, equations, and nonlinear contact
- Eigen CPU solvers with optional oneMKL/PARDISO acceleration
- Optional CUDA solvers and optional NVIDIA cuDSS
- Built-in input-deck DSL documentation

## Requirements

Always required:

- CMake 3.24 or newer
- A C++17 compiler
- Eigen, Spectra, and argparse headers

Optional components:

- Ninja for Linux/macOS presets
- Visual Studio 2022 for Windows presets
- Intel oneAPI oneMKL plus the oneAPI compiler runtime for MKL presets
- CUDA Toolkit for CUDA presets
- NVIDIA cuDSS for cuDSS presets
- GoogleTest when `FEMASTER_BUILD_TESTS=ON`

CMake first searches installed packages and then common include directories. A
custom header root containing `Eigen/`, `Spectra/`, and `argparse/` can be set
with:

```text
-DFEMASTER_DEPS_INCLUDE_DIR=/path/to/include
```

## Preset workflow

List the presets available on the current host:

```bash
cmake --list-presets
```

Configure and build with the same preset name:

```bash
cmake --preset linux-release
cmake --build --preset linux-release
```

Build presets use eight parallel jobs. CMake prints the selected compiler,
OpenMP runtime, MKL threading layer, CUDA version, cuDSS version, and runtime
linkage during configuration. Ninja and Visual Studio provide the per-file
progress output.

### Linux and WSL

CPU build with OpenMP:

```bash
cmake --preset linux-release
cmake --build --preset linux-release
```

For MKL, initialize oneAPI in the shell first:

```bash
source /opt/intel/oneapi/setvars.sh
cmake --preset linux-mkl
cmake --build --preset linux-mkl
```

`linux-mkl` uses threaded, statically linked MKL and links libstdc++/libgcc
statically. `linux-mkl-sequential` selects sequential MKL. The resulting file
is:

```text
build/linux-mkl/bin/FEMaster
```

CUDA variants are:

```bash
cmake --preset linux-cuda
cmake --build --preset linux-cuda

cmake --preset linux-cuda-mkl
cmake --build --preset linux-cuda-mkl
```

For cuDSS, set its installation root before configuring:

```bash
export CUDSS_DIR=/path/to/libcudss-linux-x86_64
cmake --preset linux-cuda-cudss
cmake --build --preset linux-cuda-cudss
```

`linux-full` enables CUDA, cuDSS, and oneMKL together.

### Windows

Run CMake from an **Intel oneAPI command prompt for Visual Studio 2022** when
using MKL. A normal Visual Studio x64 prompt is sufficient without MKL.

```powershell
cmake --preset windows-mkl
cmake --build --preset windows-mkl
```

The executable and license files are written to:

```text
build\windows-mkl\bin\Release\FEMaster.exe
```

The threaded static-MKL preset still needs Intel's OpenMP runtime DLL. CMake
finds `libiomp5md.dll` in the oneAPI installation and copies it beside the
executable automatically. Use `windows-mkl-sequential` if that runtime is not
wanted.

CUDA variants are configured from a prompt where the CUDA Toolkit is visible:

```powershell
cmake --preset windows-cuda
cmake --build --preset windows-cuda

cmake --preset windows-cuda-mkl
cmake --build --preset windows-cuda-mkl
```

Set `CUDSS_DIR` to the unpacked cuDSS root for `windows-cuda-cudss`.
`windows-full` enables CUDA, cuDSS, and oneMKL together.

### macOS

macOS provides CPU builds:

```bash
cmake --preset macos-release
cmake --build --preset macos-release
```

After installing an OpenMP runtime discoverable by CMake, use
`macos-openmp`. Current NVIDIA CUDA toolkits do not support macOS, so CMake
rejects CUDA on that platform intentionally. No macOS MKL preset is provided;
the supported MKL targets for this project are Linux and Windows.

## CMake options

| Option | Default | Meaning |
|---|---:|---|
| `FEMASTER_ENABLE_OPENMP` | `OFF` | Enable FEMaster OpenMP loops |
| `FEMASTER_ENABLE_MKL` | `OFF` | Enable oneMKL/PARDISO |
| `FEMASTER_MKL_SEQUENTIAL` | `OFF` | Use sequential instead of threaded MKL |
| `FEMASTER_MKL_STATIC` | `ON` | Link MKL component libraries statically |
| `FEMASTER_ENABLE_CUDA` | `OFF` | Build CUDA solver backends |
| `FEMASTER_ENABLE_CUDSS` | `OFF` | Use cuDSS; requires CUDA |
| `FEMASTER_DOUBLE_PRECISION` | `ON` | Use double precision on CPU and CUDA |
| `FEMASTER_STATIC_RUNTIME` | `OFF` | Statically link compiler C/C++ runtimes |
| `FEMASTER_FULLY_STATIC` | `OFF` | Request a fully static Linux executable |
| `FEMASTER_BUILD_TESTS` | `OFF` | Build GoogleTest tests |
| `FEMASTER_TIME_REPORT` | `OFF` | Enable compiler timing diagnostics |

Options can override a preset:

```bash
cmake --preset linux-release -DFEMASTER_DOUBLE_PRECISION=OFF
cmake --build --preset linux-release
```

For CUDA architecture selection, use CMake's standard setting, for example:

```text
-DCMAKE_CUDA_ARCHITECTURES=86
```

`FEMASTER_FULLY_STATIC` is separate from the distributable MKL presets. Fully
static glibc builds depend on static versions of every system library and are
incompatible with CUDA.

## Runtime distribution

CPU presets without dynamic third-party libraries need no FEMaster-specific
runtime installation. The static MKL presets embed MKL itself. On Windows,
threaded MKL additionally ships the automatically copied `libiomp5md.dll`.

CUDA cannot be folded into a single portable executable: CUDA builds require a
compatible NVIDIA driver and may require CUDA/cuDSS runtime libraries according
to NVIDIA's redistribution rules. Windows cuDSS builds copy the imported cuDSS
DLL beside FEMaster; other NVIDIA runtime dependencies remain external.

## Running FEMaster

```bash
build/linux-release/bin/FEMaster model.inp --ncpus 8
```

```powershell
.\build\windows-mkl\bin\Release\FEMaster.exe .\model.inp --ncpus 8
```

The input filename is positional. Results default to the same path with a
`.res` suffix; `--output` selects another result path.

DSL documentation is available from the executable:

```bash
build/linux-release/bin/FEMaster --document --doc-list
build/linux-release/bin/FEMaster --document --doc-show LOADCASE
build/linux-release/bin/FEMaster --document --doc-search shell
```

## Tests

Enable tests on any suitable preset, then build and run CTest:

```bash
cmake --preset linux-release -DFEMASTER_BUILD_TESTS=ON
cmake --build --preset linux-release
ctest --test-dir build/linux-release --output-on-failure
```

For Visual Studio builds, add `-C Release` to the CTest command.

## License

FEMaster is licensed under the MIT License; see [LICENSE.txt](LICENSE.txt).
Third-party components retain their own licenses; see
[THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt) and `licenses/`.
