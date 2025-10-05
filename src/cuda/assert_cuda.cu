/******************************************************************************
 * @file assert_cuda.cu
 * @brief Implements CUDA error translation helpers.
 *
 * @see src/cuda/assert_cuda.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#ifdef SUPPORT_GPU

#include "cuda.h"

#include <cstdio>
#include <iostream>
#include <string>

namespace fem {
namespace cuda {

void gpu_assert(cudaError_t code, const char* file, int line, bool abort) {
    if (code != cudaSuccess) {
        std::fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) {
            std::exit(static_cast<int>(code));
        }
    }
}

void gpu_assert(cusolverStatus_t code, const char* file, int line, bool abort) {
    if (code != CUSOLVER_STATUS_SUCCESS) {
        std::string label;
        switch (code) {
            case CUSOLVER_STATUS_SUCCESS: label = "CUSOLVER_STATUS_SUCCESS"; break;
            case CUSOLVER_STATUS_NOT_INITIALIZED: label = "CUSOLVER_STATUS_NOT_INITIALIZED"; break;
            case CUSOLVER_STATUS_ALLOC_FAILED: label = "CUSOLVER_STATUS_ALLOC_FAILED"; break;
            case CUSOLVER_STATUS_INVALID_VALUE: label = "CUSOLVER_STATUS_INVALID_VALUE"; break;
            case CUSOLVER_STATUS_ARCH_MISMATCH: label = "CUSOLVER_STATUS_ARCH_MISMATCH"; break;
            case CUSOLVER_STATUS_EXECUTION_FAILED: label = "CUSOLVER_STATUS_EXECUTION_FAILED"; break;
            case CUSOLVER_STATUS_INTERNAL_ERROR: label = "CUSOLVER_STATUS_INTERNAL_ERROR"; break;
            case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
                label = "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
                break;
            default: label = "<unknown>"; break;
        }
        std::fprintf(stderr, "GPUassert: %s %s %d\n", label.c_str(), file, line);
        if (abort) {
            std::exit(static_cast<int>(code));
        }
    }
}

void gpu_assert(cublasStatus_t code, const char* file, int line, bool abort) {
    if (code != CUBLAS_STATUS_SUCCESS) {
        std::string label;
        switch (code) {
            case CUBLAS_STATUS_SUCCESS: label = "CUBLAS_STATUS_SUCCESS"; break;
            case CUBLAS_STATUS_NOT_INITIALIZED: label = "CUBLAS_STATUS_NOT_INITIALIZED"; break;
            case CUBLAS_STATUS_ALLOC_FAILED: label = "CUBLAS_STATUS_ALLOC_FAILED"; break;
            case CUBLAS_STATUS_INVALID_VALUE: label = "CUBLAS_STATUS_INVALID_VALUE"; break;
            case CUBLAS_STATUS_ARCH_MISMATCH: label = "CUBLAS_STATUS_ARCH_MISMATCH"; break;
            case CUBLAS_STATUS_MAPPING_ERROR: label = "CUBLAS_STATUS_MAPPING_ERROR"; break;
            case CUBLAS_STATUS_EXECUTION_FAILED: label = "CUBLAS_STATUS_EXECUTION_FAILED"; break;
            case CUBLAS_STATUS_INTERNAL_ERROR: label = "CUBLAS_STATUS_INTERNAL_ERROR"; break;
            case CUBLAS_STATUS_NOT_SUPPORTED: label = "CUBLAS_STATUS_NOT_SUPPORTED"; break;
            default: label = "<unknown>"; break;
        }
        std::fprintf(stderr, "GPUassert: %s %s %d\n", label.c_str(), file, line);
        if (abort) {
            std::exit(static_cast<int>(code));
        }
    }
}

void gpu_assert(cusparseStatus_t code, const char* file, int line, bool abort) {
    if (code != CUSPARSE_STATUS_SUCCESS) {
        std::fprintf(stderr, "GPUassert: %s %s %d\n", cusparseGetErrorString(code), file, line);
        if (abort) {
            std::exit(static_cast<int>(code));
        }
    }
}

} // namespace cuda
} // namespace fem

#endif
