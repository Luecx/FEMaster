/**
 * @file cuda_vec.h
 * @brief Binds owned device-vector storage to a cuSPARSE dense-vector descriptor.
 *
 * CudaArray owns allocation, transfers and release of the device buffer.
 * CudaVector adds the descriptor required by sparse GPU matrix-vector products
 * and exposes it to solver backends. The descriptor and scalar storage use the
 * configured CudaPrecision throughout their common lifetime.
 *
 * @see CudaArray
 * @see Manager
 */

#pragma once

#include "../core/types_num.h"

#include "assert_cuda.h"
#include "cuda.h"
#include "cuda_array.h"
#include "cuda_defs.h"

#ifdef SUPPORT_GPU
namespace fem::cuda {
/**
 * @brief Owns a device vector and its matching cuSPARSE view.
 *
 * Construction allocates the inherited CudaArray buffer and creates a dense
 * descriptor of the same length and precision. Destruction releases the
 * descriptor before the base releases its storage; cuSPARSE does not own the
 * underlying allocation. Solver kernels borrow the descriptor through the
 * conversion operator while this object and its buffer remain alive.
 */
struct CudaVector : CudaArray<CudaPrecision> {
    cusparseDnVecDescr_t descr;

    explicit CudaVector(int size)
        : CudaArray(size) {
        runtime_check_cuda(cusparseCreateDnVec(&descr , size, this->data , CUDA_P_TYPE));
    }

    ~CudaVector() {
        runtime_check_cuda(cusparseDestroyDnVec(descr));
    }

    operator cusparseDnVecDescr_t() {
        return descr;
    }
};
}    // namespace cuda

#endif
