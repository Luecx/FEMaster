
#ifdef SUPPORT_GPU

#include "cuda.h"

#include <iostream>
#include <string>

namespace cuda {

Manager manager;

void Manager::create_cuda() {
    if (loaded)
        return;

    int device_count;
    cudaGetDeviceCount(&device_count);

    // check that there is atleast 1 device
    if(device_count == 0){
        log_warning(device_count, "No devices found");
        return;
    }

    // loading handles
    cusolverSpCreate(&handle_cusolver_sp);
    cusolverDnCreate(&handle_cusolver_dn);
    cublasCreate_v2 (&handle_cublas);
    cusparseCreate  (&handle_cusparse);

    // get properties
    cudaGetDeviceProperties(&gpu_properties, 0);
    cudaSetDevice(0);

    // set loading flag
    loaded = true;

    int major=-1,minor=-1,patch=-1;
    cusolverGetProperty(MAJOR_VERSION, &major);
    cusolverGetProperty(MINOR_VERSION, &minor);
    cusolverGetProperty(PATCH_LEVEL  , &patch);

    // logging
    log_info(loaded, "Creating CUDA handles");
    log_info(loaded, "Using " + std::string(gpu_properties.name));
    log_info(loaded, "Memory: " + std::to_string(gpu_properties.totalGlobalMem));
    log_info(loaded, "CuSolver: " + std::to_string(major) + "." + std::to_string(minor) + "." + std::to_string(patch));
}

void Manager::destroy_cuda() {
    if (!(loaded))
        return;
    cusolverSpDestroy(handle_cusolver_sp);
    cusolverDnDestroy(handle_cusolver_dn);
    cublasDestroy_v2 (handle_cublas);
    cusparseDestroy  (handle_cusparse);
    loaded = false;
}

bool Manager::is_loaded() {
    return loaded;
}

Manager::~Manager() {
    destroy_cuda();
}

size_t Manager::mem_free() {
    size_t m_free, m_total;
    runtime_check_cuda(cudaMemGetInfo(&m_free, &m_total));
    return m_free;
}

}    // namespace cuda
#endif