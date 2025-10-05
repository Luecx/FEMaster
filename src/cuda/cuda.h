/******************************************************************************
 * @file cuda.h
 * @brief Declares the CUDA resource manager and associated helpers.
 *
 * Encapsulates handle creation and destruction for the cuSOLVER, cuBLAS,
 * cuSPARSE, and cuRAND libraries.
 *
 * @see src/cuda/cuda.cu
 * @see src/cuda/assert_cuda.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#ifdef SUPPORT_GPU

#include "../core/core.h"

#include <cublas_v2.h>
#include <cuComplex.h>
#include <curand.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>

namespace fem {
namespace cuda {

/******************************************************************************
 * @struct Manager
 * @brief Owns CUDA library handles and exposes GPU properties.
 ******************************************************************************/
struct Manager {
    cusolverSpHandle_t handle_cusolver_sp{};
    cusolverDnHandle_t handle_cusolver_dn{};
    cublasHandle_t handle_cublas{};
    cusparseHandle_t handle_cusparse{};

    ~Manager();

    void create_cuda();
    void destroy_cuda();

    std::size_t mem_free();
    bool is_loaded();

private:
    bool loaded = false;
    cudaDeviceProp gpu_properties{};
};

extern Manager manager;

} // namespace cuda
} // namespace fem

#endif
