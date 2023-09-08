#pragma once

#ifdef SUPPORT_GPU

#include "../core/core.h"

#include <cuComplex.h>
#include <cublas_v2.h>
#include <curand.h>
#include <cusolverDn.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>

namespace cuda {

struct Manager {

    // handles
    cusolverSpHandle_t handle_cusolver_sp;
    cusolverDnHandle_t handle_cusolver_dn;
    cublasHandle_t     handle_cublas;
    cusparseHandle_t   handle_cusparse;

    private:
    // state of cuda (loaded or not)
    bool loaded = false;

    // data about graphics card, extracted on loading
    cudaDeviceProp gpu_properties;

    public:
    void   destroy_cuda();
    void   create_cuda();

    size_t mem_free();
    bool   is_loaded();

    ~Manager();
};

extern Manager manager;

}    // namespace cuda

#endif