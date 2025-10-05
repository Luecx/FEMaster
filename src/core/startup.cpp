/******************************************************************************
 * @file startup.cpp
 * @brief Implements the startup helper that reports build information.
 *
 * @see src/core/startup.h
 * @see src/core/version.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "startup.h"

#include "core.h"
#include "version.h"

#include <Spectra/Util/Version.h>
#include <iostream>
#include <string>

namespace fem {
namespace startup {
namespace {

void print_banner() {
    std::cout << "**********************************************************************\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*                         FEMaster                                   *\n";
    std::cout << "*                          v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH <<
        "                                    *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Copyright (c) 2024, Finn Eggers                          *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           All rights reserved.                                     *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           This program is provided \"as is\" without any             *\n";
    std::cout << "*           warranty of any kind, either expressed or                *\n";
    std::cout << "*           implied. Use at your own risk.                           *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Decompilation or reverse engineering of this             *\n";
    std::cout << "*           software is strictly prohibited.                         *\n";
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
    std::cout << "*           Bytes of GPU Precision: N/A                             *\n";
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
    std::cout << "*                                                                    *\n";
    std::cout << "*           FEMaster Version      : " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH <<
        "                            *\n";
    std::cout << "*           Eigen Version         : " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." <<
        EIGEN_MINOR_VERSION << "                           *\n";
    std::cout << "*           Spectra Version       : " << SPECTRA_MAJOR_VERSION << "." << SPECTRA_MINOR_VERSION << "." <<
        SPECTRA_PATCH_VERSION << "                            *\n";
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
