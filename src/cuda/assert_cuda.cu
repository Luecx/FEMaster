
#ifdef SUPPORT_GPU

#include "cuda.h"

#include <cstdio>
#include <iostream>

namespace cuda {

void gpu_assert(cudaError_t code, const char* file, int line, bool abort) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
}

void gpu_assert(cusolverStatus_t code, const char* file, int line, bool abort) {
    if (code != CUSOLVER_STATUS_SUCCESS) {
        std::string s_code;
        switch (code) {
            case CUSOLVER_STATUS_SUCCESS: s_code = "CUSOLVER_STATUS_SUCCESS"; break;
            case CUSOLVER_STATUS_NOT_INITIALIZED: s_code = "CUSOLVER_STATUS_NOT_INITIALIZED"; break;
            case CUSOLVER_STATUS_ALLOC_FAILED: s_code = "CUSOLVER_STATUS_ALLOC_FAILED"; break;
            case CUSOLVER_STATUS_INVALID_VALUE: s_code = "CUSOLVER_STATUS_INVALID_VALUE"; break;
            case CUSOLVER_STATUS_ARCH_MISMATCH: s_code = "CUSOLVER_STATUS_ARCH_MISMATCH"; break;
            case CUSOLVER_STATUS_EXECUTION_FAILED: s_code = "CUSOLVER_STATUS_EXECUTION_FAILED"; break;
            case CUSOLVER_STATUS_INTERNAL_ERROR: s_code = "CUSOLVER_STATUS_INTERNAL_ERROR"; break;
            case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED: s_code = "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED"; break;
            default: s_code = "<unknown>";
        }

        fprintf(stderr, "GPUassert: %s %s %d\n", s_code.c_str(), file, line);
        if (abort)
            exit(code);
    }
}

void gpu_assert(cublasStatus_t code, const char* file, int line, bool abort) {
    if (code != CUBLAS_STATUS_SUCCESS) {
        std::string s_code;
        switch (code) {
            case CUBLAS_STATUS_SUCCESS: s_code = "CUBLAS_STATUS_SUCCESS"; break;
            case CUBLAS_STATUS_NOT_INITIALIZED: s_code = "CUBLAS_STATUS_NOT_INITIALIZED"; break;
            case CUBLAS_STATUS_ALLOC_FAILED: s_code = "CUBLAS_STATUS_ALLOC_FAILED"; break;
            case CUBLAS_STATUS_INVALID_VALUE: s_code = "CUBLAS_STATUS_INVALID_VALUE"; break;
            case CUBLAS_STATUS_ARCH_MISMATCH: s_code = "CUBLAS_STATUS_ARCH_MISMATCH"; break;
            case CUBLAS_STATUS_MAPPING_ERROR: s_code = "CUBLAS_STATUS_MAPPING_ERROR"; break;
            case CUBLAS_STATUS_EXECUTION_FAILED: s_code = "CUBLAS_STATUS_EXECUTION_FAILED"; break;
            case CUBLAS_STATUS_INTERNAL_ERROR: s_code = "CUBLAS_STATUS_INTERNAL_ERROR"; break;
            case CUBLAS_STATUS_NOT_SUPPORTED: s_code = "CUBLAS_STATUS_NOT_SUPPORTED"; break;
            default: s_code = "<unknown>";
        }
        fprintf(stderr, "GPUassert: %s %s %d\n", s_code.c_str(), file, line);
        if (abort)
            exit(code);
    }
}

void gpu_assert(cusparseStatus_t code, const char* file, int line, bool abort) {
    if (code != CUSPARSE_STATUS_SUCCESS) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cusparseGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
}
}    // namespace cuda
#endif