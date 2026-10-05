/**
 * @file startup.cpp
 * @brief Implements the startup helper that reports build information.
 *
 * The core startup object prints solver and dependency versions together with
 * compile-time backend availability before main. Version values occupy fixed
 * output columns so their digit counts do not shift the banner border.
 *
 * @see src/core/startup.h
 * @see src/core/version.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "startup.h"

#include "core.h"
#include "version.h"

#include <Spectra/Util/Version.h>
#include <iomanip>
#include <iostream>
#include <string>

namespace fem {
namespace startup {
namespace {

/**
 * Prints the solver banner and compile-time build information.
 *
 * Version strings are padded on the right to keep the closing border in column
 * 70 regardless of their digit counts within the available width. Each version
 * field restores right alignment afterward to preserve subsequent formatting.
 * The global Startup instance invokes this output before command-line parsing.
 */
void print_banner() {
    // Assemble complete version values before applying output field widths.
    const std::string solver_version  = std::to_string(VERSION_MAJOR) + "." +
                                        std::to_string(VERSION_MINOR) + "." +
                                        std::to_string(VERSION_PATCH);
    const std::string eigen_version   = std::to_string(EIGEN_WORLD_VERSION) + "." +
                                        std::to_string(EIGEN_MAJOR_VERSION) + "." +
                                        std::to_string(EIGEN_MINOR_VERSION);
    const std::string spectra_version = std::to_string(SPECTRA_MAJOR_VERSION) + "." +
                                        std::to_string(SPECTRA_MINOR_VERSION) + "." +
                                        std::to_string(SPECTRA_PATCH_VERSION);

    // Reserve the remaining columns for each value and its trailing padding.
    std::cout << "**********************************************************************\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*                         FEMaster                                   *\n";
    std::cout << "*                          v" << std::left << std::setw(41)
              << solver_version << std::right << "*\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Copyright (c) 2024, Finn Eggers                            *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Licensed under the MIT License.                          *\n";
    std::cout << "*           See LICENSE.txt for the complete license terms.          *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           This program is provided \"as is\" without any             *\n";
    std::cout << "*           warranty of any kind, either expressed or                *\n";
    std::cout << "*           implied. Use at your own risk.                           *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*....................................................................*\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Build Information                                        *\n";
    std::cout << "*           Bytes of CPU Precision: " << sizeof(Precision);
    std::cout << std::string(33 - std::to_string(sizeof(Precision)).length(), ' ') << "*\n";
#ifdef SUPPORT_GPU
    std::cout << "*           Bytes of GPU Precision: " << sizeof(CudaPrecision);
    std::cout << std::string(33 - std::to_string(sizeof(CudaPrecision)).length(), ' ') << "*\n";
#else
    std::cout << "*           Bytes of GPU Precision: N/A                              *\n";
#endif
#ifdef SUPPORT_GPU
    std::cout << "*           GPU Supported         : Yes                              *\n";
#else
    std::cout << "*           GPU Supported         : No                               *\n";
#endif
#ifdef _OPENMP
    std::cout << "*           OPENMP Supported      : Yes                              *\n";
#else
    std::cout << "*           OPENMP Supported      : No                               *\n";
#endif
#ifdef USE_MKL
    std::cout << "*           MKL Supported         : Yes                              *\n";
#else
    std::cout << "*           MKL Supported         : No                               *\n";
#endif
#ifdef USE_CUDSS
    std::cout << "*           cuDSS Supported       : Yes                              *\n";
#else
    std::cout << "*           cuDSS Supported       : No                               *\n";
#endif
    std::cout << "*                                                                    *\n";
    std::cout << "*           FEMaster Version      : " << std::left << std::setw(33)
              << solver_version << std::right << "*\n";
    std::cout << "*           Eigen Version         : " << std::left << std::setw(33)
              << eigen_version << std::right << "*\n";
    std::cout << "*           Spectra Version       : " << std::left << std::setw(33)
              << spectra_version << std::right << "*\n";
    std::cout << "*                                                                    *\n";
    std::cout << "**********************************************************************\n";
}
} // namespace

Startup::Startup() {
    print_banner();
}

Startup instance{};
} // namespace startup
} // namespace fem
