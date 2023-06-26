#pragma once

#include "assert_cuda.h"
#include "cuda.h"
#include "cuda_array.h"
#include "cuda_defs.h"

#ifdef SUPPORT_GPU
namespace cuda {

struct CudaVector : CudaArray<Precision> {
    cusparseDnVecDescr_t descr;

    explicit CudaVector(int size)
        : CudaArray(size) {
        runtime_check_cuda(cusparseCreateDnVec(&descr , size, this->data , CUDA_P_TYPE));
    }

    ~CudaVector() override {
        runtime_check_cuda(cusparseDestroyDnVec(descr));
    }

    operator cusparseDnVecDescr_t(){
        return descr;
    }
};
}    // namespace cuda

#endif
