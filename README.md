<div align="center">

# FEMaster

### Finite Elements in Mechanical Analysis and Simulation for Technical Engineering Research

Open-source finite-element solver developed by [Finn Eggers](https://github.com/Luecx).

[![Release](https://img.shields.io/github/v/release/Luecx/FEMaster)](https://github.com/Luecx/FEMaster/releases)

[Wiki](https://github.com/Luecx/FEMaster/wiki) ·
[Documentation](https://github.com/Luecx/FEMaster/tree/master/documentation) ·
[Releases](https://github.com/Luecx/FEMaster/releases) ·
[Contributors](https://github.com/Luecx/FEMaster/graphs/contributors) ·
[Issues](https://github.com/Luecx/FEMaster/issues)

</div>

---

FEMaster is a modern C++ finite-element solver for structural mechanics.
It supports linear and nonlinear analyses across Linux, Windows, and macOS.

## Features

- Linear and geometrically nonlinear static analysis
- Eigenfrequency and linear buckling analysis
- Implicit transient analysis
- Beam, shell, solid, truss, and point-mass elements
- Supports, couplings, equations, ties, connectors, and nonlinear contact
- Abaqus-style input decks
- Parallel CPU execution with OpenMP
- Intel oneMKL / PARDISO acceleration
- NVIDIA CUDA and cuDSS solver backends
- Apple Accelerate backend on macOS

## Getting started

Compilation instructions, platform-specific requirements, input syntax, examples,
and further documentation are available in the
[FEMaster Wiki](https://github.com/Luecx/FEMaster/wiki).

The full FEMaster user manual is available in the
[documentation](https://github.com/Luecx/FEMaster/tree/master/documentation).

## License

FEMaster is distributed under the [MIT License](LICENSE.txt).

Third-party components retain their respective licenses. See
[THIRD_PARTY_NOTICES.txt](THIRD_PARTY_NOTICES.txt) and the `licenses/` directory
for details.

