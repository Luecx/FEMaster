#pragma once

#ifdef SUPPORT_GPU

#include <cuComplex.h>
#include <cublas_v2.h>
#include <curand.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>
#include "cuda.h"

#include <iostream>

namespace cuda{

void gpu_assert(cudaError_t      code, const char* file, int line, bool abort = true);
void gpu_assert(cusolverStatus_t code, const char* file, int line, bool abort = true);
void gpu_assert(cusparseStatus_t code, const char* file, int line, bool abort = true);
void gpu_assert(cublasStatus_t   code, const char* file, int line, bool abort = true);

}
#define runtime_check_cuda(ans)                                                                                        \
    { cuda::gpu_assert((ans), __FILE__, __LINE__); }

#endif