//
// Created by Luecx on 06.09.2023.
//

#include "startup.h"
#include "core.h"
#include <iostream>


void copyright() {
    std::cout << "**********************************************************************\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*                         FEMaster                                   *\n";
    std::cout << "*                                                                    *\n";
    std::cout << "*           Copyright (c) 2023, Finn Eggers                          *\n";
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
    std::cout << "*           Bytes of Precision    : " << sizeof(Precision);
    std::cout << std::string(33 - std::to_string(sizeof(Precision)).length(), ' ') << "*\n";
    std::cout << "*           Bytes of CudaPrecision: " << sizeof(CudaPrecision);
    std::cout << std::string(33 - std::to_string(sizeof(CudaPrecision)).length(), ' ') << "*\n";
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
    std::cout << "**********************************************************************\n";
}

Startup::Startup() {
    copyright();
}

Startup startup{};
