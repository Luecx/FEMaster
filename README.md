# FEMaster

**FEMaster** is a powerful tool designed for solving linear finite element (FE) problems. It supports linear analysis and linear frequency extraction using high-performance algorithms and optional GPU acceleration. FEMaster offers various customization features, including support for multi-core processing via OpenMP, optional integration with the Math Kernel Library (MKL) for optimized performance, and CUDA for leveraging GPU computations. This tool is ideal for engineers and researchers looking to perform fast and efficient FE analysis on large-scale problems.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Optional Features](#optional-features)
- [Usage](#usage)
    - [Compiling the Project](#compiling-the-project)
    - [Compilation Options](#compilation-options)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)

## Features
- **Linear Analysis**: Solve linear finite element problems efficiently.
- **Linear Frequency Extraction**: Perform frequency extraction for linear systems.
- **CUDA Support**: Optional GPU acceleration with CUDA.
- **OpenMP**: Enable multi-threaded execution for faster computations.
- **MKL Integration**: Use Intel's Math Kernel Library (MKL) for optimized performance.
- **Cross-platform Support**: Compatible with Linux, macOS, and Windows.

## Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/Luecx/FEMaster.git
cd FEMaster
```

### Step 2: Install Required Dependencies
This project depends on several libraries such as Eigen (for matrix operations) and Spectra (for eigenvalue problems). The following script will download and install these libraries.

```bash
#!/bin/bash

# Create a directory to store libraries
mkdir -p ~/libs

# Install Eigen
cd ~/libs
git clone https://gitlab.com/libeigen/eigen.git
sudo cp -r eigen/Eigen /usr/local/include/

# Install Spectra
git clone https://github.com/yixuan/spectra.git
sudo cp -r spectra/include/* /usr/local/include/
```

Save the above script as `install_dependencies.sh`, make it executable, and run it:

```bash
chmod +x install_dependencies.sh
./install_dependencies.sh
```

### Step 3: (Optional) Install MKL for Optimal Performance
Intel's Math Kernel Library (MKL) is optional but provides optimized performance for numerical computations. To install MKL, follow the instructions on the official [Intel MKL website](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html).

Once installed, set the `MKLROOT` environment variable:

```bash
source /opt/intel/oneapi/setvars.sh
```

### Step 4: (Optional) Install CUDA
For GPU acceleration, you can install the CUDA toolkit. Download and install it from the official [NVIDIA CUDA website](https://developer.nvidia.com/cuda-toolkit).

## Dependencies

FEMaster relies on the following libraries for core functionality:
- **Eigen**: For matrix operations (download from [Eigen website](https://eigen.tuxfamily.org/dox/))
- **Spectra**: For eigenvalue computations (download from [Spectra GitHub](https://github.com/yixuan/spectra))

Optional dependencies:
- **MKL**: Intel's Math Kernel Library for optimized performance.
- **CUDA**: NVIDIA's toolkit for GPU acceleration.

## Optional Features

- **OpenMP**: Enabled by default for multi-core processing.
- **MKL**: Can be optionally enabled for high-performance linear algebra routines.
- **CUDA Double Precision**: Can be enabled for high-precision GPU computation.

## Usage

### Compiling the Project

FEMaster can be compiled in different modes depending on the desired features (CPU-only, GPU acceleration, etc.). You can customize the build with flags for OpenMP, MKL, and CUDA support.

#### General Compilation Command
```bash
make all
```

This command will compile both the CPU-only and GPU-accelerated versions of the code.

#### Building CPU-Only Version
```bash
make cpu
```

#### Building GPU Version (with CUDA support)
```bash
make gpu
```

#### Running the Tests
```bash
make tests
```

#### Clean Compiled Files
```bash
make clean
```

### Compilation Options

Several build flags allow you to enable or disable features as needed:

| Flag        | Default | Description                                              |
|-------------|---------|----------------------------------------------------------|
| `openmp`    | 1       | Enables/disables OpenMP for multi-core support            |
| `mkl`       | 0       | Enables/disables MKL integration for optimized performance |
| `cuda_dp`   | 1       | Enables/disables CUDA double precision                    |
| `ar_pcs`    | 0       | Enables/disables array process debug information          |
| `debug`     | 0       | Enables/disables debug mode                               |

You can override these options during compilation by passing them as environment variables. For example, to disable OpenMP and enable MKL:

```bash
make cpu openmp=0 mkl=1
```

### Sample Compilation Process
```bash
make cpu            # Compile for CPU execution
make gpu            # Compile for GPU execution (CUDA required)
make all mkl=1      # Compile with MKL support enabled
```

## Testing

To run the test suite, use the following command:
```bash
make tests
```

The tests are built using Google Test, and they ensure that the functionality of the FEMaster code is working correctly across different configurations.

## Contributing

Contributions are welcome! Please submit issues and pull requests to the repository, and feel free to suggest improvements or new features.

## License

FEMaster is licensed under the MIT License. See `LICENSE` for more details.