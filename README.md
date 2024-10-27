# FEMaster [![CodeFactor](https://www.codefactor.io/repository/github/luecx/femaster/badge)](https://www.codefactor.io/repository/github/luecx/femaster)

**FEMaster** is a powerful tool designed for solving linear finite element (FE) problems. It supports linear analysis and linear frequency extraction using high-performance algorithms and optional GPU acceleration. FEMaster offers various customization features, including support for multi-core processing via OpenMP, optional integration with the Math Kernel Library (MKL) for optimized performance, and CUDA for leveraging GPU computations. This tool is ideal for engineers and researchers looking to perform fast and efficient FE analysis on large-scale problems.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Optional Features](#optional-features)
- [Usage](#usage)
  - [Compiling the Project](#compiling-the-project)
  - [Build Configuration](#build-configuration)
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
- **Cross-platform Support**: Compatible with Linux, macOS, and Windows (including Apple Silicon support).

## Installation

### Step 1: Clone the Repository
```bash
git clone https://github.com/Luecx/FEMaster.git
cd FEMaster
```

### Step 2: Install Required Dependencies
This project depends on several libraries such as Eigen (for matrix operations) and Spectra (for eigenvalue problems). It also makes use of argparse for command line arguments. The following script will download and install these libraries.

```bash
#!/bin/bash

# Create a directory to store libraries
sudo mkdir -p ~/libs
sudo mkdir -p /usr/local/include

# Install Eigen
cd ~/libs
git clone https://gitlab.com/libeigen/eigen.git
sudo cp -r eigen/Eigen /usr/local/include/

# Install Spectra
git clone https://github.com/yixuan/spectra.git
sudo cp -r spectra/include/Spectra /usr/local/include/

# install argparse
git clone https://github.com/p-ranav/argparse.git
sudo cp -r argparse/include/argparse /usr/local/include/
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

### Required Dependencies
- **Eigen**: For matrix operations (download from [Eigen website](https://eigen.tuxfamily.org/dox/))
- **Spectra**: For eigenvalue computations (download from [Spectra GitHub](https://github.com/yixuan/spectra))

### Optional Dependencies
- **MKL**: Intel's Math Kernel Library for optimized performance
- **CUDA**: NVIDIA's toolkit for GPU acceleration
- **Google Test**: For running the test suite

## Optional Features

The project supports several optional features that can be enabled or disabled during compilation:

| Feature | Description |
|---------|-------------|
| OpenMP | Multi-threading support for parallel computation |
| MKL | Intel Math Kernel Library integration |
| Sequential MKL | Single-threaded MKL operation |
| CUDA Double Precision | High-precision GPU computation |
| Double Precision | Double-precision floating-point computation |
| Debug Mode | Enhanced debugging information |

## Usage

### Compiling the Project

FEMaster provides several make targets for different build configurations:

```bash
make all      # Build everything (CPU, GPU, and tests)
make cpu      # Build CPU-only version
make gpu      # Build GPU-enabled version
make tests    # Build test suite
make clean    # Remove build artifacts
make info     # Display build configuration
make help     # Show available make targets
```

### Build Configuration

You can view your current build configuration using:

```bash
make info
```

This will display information about:
- Project name
- Platform and architecture
- Compiler version
- Compiler flags
- Feature flags (MKL, OpenMP, Debug mode)

### Compilation Options

The following build flags can be configured during compilation:

| Flag | Default | Description |
|------|---------|-------------|
| `openmp` | 1 | Enable/disable OpenMP support |
| `mkl` | 0 | Enable/disable MKL integration |
| `mkl_sequential` | 0 | Enable/disable Sequential MKL |
| `cuda_dp` | 1 | Enable/disable CUDA double precision |
| `debug` | 0 | Enable/disable debug mode |
| `double_precision` | 1 | Enable/disable double precision |

To use these flags, append them to your make command:

```bash
make cpu openmp=0 mkl=1              # Build CPU version without OpenMP but with MKL
make gpu cuda_dp=1 double_precision=1 # Build GPU version with double precision
make all debug=1                     # Build everything in debug mode
```

### Platform-Specific Features

The build system automatically detects and adapts to different platforms:
- Supports both x86_64 and ARM64 architectures
- Automatic detection of Clang vs GCC
- Platform-specific library paths (especially for Apple Silicon)
- Appropriate MKL configuration based on platform

## Testing

The project uses Google Test for its test suite. To run the tests:

```bash
make tests
./bin/FEMaster_test
```

Test binaries are automatically configured for your platform, including proper library paths for both Intel and Apple Silicon machines.

## Contributing

Contributions are welcome! Please submit issues and pull requests to the repository, and feel free to suggest improvements or new features.

## License

FEMaster is licensed under the MIT License. See `LICENSE` for more details.