/******************************************************************************
 * @file assert_cuda.h
 * @brief Declares helper functions that translate CUDA error codes into exceptions.
 *
 * The helpers provide consistent diagnostics for CUDA runtime, cuSOLVER,
 * cuSPARSE, and cuBLAS calls, integrating with the CPU-side error handling.
 *
 * @see src/cuda/assert_cuda.cu
 * @see src/cuda/cuda.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#ifdef SUPPORT_GPU

#include <cublas_v2.h>
#include <cuComplex.h>
#include <curand.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>

#include <iostream>

#include "cuda.h"

namespace fem {
namespace cuda {

void gpu_assert(cudaError_t code, const char* file, int line, bool abort = true);
void gpu_assert(cusolverStatus_t code, const char* file, int line, bool abort = true);
void gpu_assert(cusparseStatus_t code, const char* file, int line, bool abort = true);
void gpu_assert(cublasStatus_t code, const char* file, int line, bool abort = true);

} // namespace cuda
} // namespace fem

#define runtime_check_cuda(ans)                                                                                        \
    { fem::cuda::gpu_assert((ans), __FILE__, __LINE__); }

#endif
