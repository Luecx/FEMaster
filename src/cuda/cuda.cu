/******************************************************************************
 * @file cuda.cu
 * @brief Implements the CUDA resource manager.
 *
 * @see src/cuda/cuda.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#ifdef SUPPORT_GPU

#include "cuda.h"
#include "assert_cuda.h"

#include <iostream>
#include <string>

namespace fem {
namespace cuda {

Manager manager;

void Manager::create_cuda() {
    if (loaded) {
        return;
    }

    int device_count = 0;
    runtime_check_cuda(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        logging::warning(false, "No CUDA devices detected");
        return;
    }

    runtime_check_cuda(cudaSetDevice(0));
    runtime_check_cuda(cudaGetDeviceProperties(&gpu_properties, 0));

    runtime_check_cuda(cusolverSpCreate(&handle_cusolver_sp));
    runtime_check_cuda(cusolverDnCreate(&handle_cusolver_dn));
    runtime_check_cuda(cublasCreate_v2(&handle_cublas));
    runtime_check_cuda(cusparseCreate(&handle_cusparse));

    int major = -1;
    int minor = -1;
    int patch = -1;
    cusolverGetProperty(MAJOR_VERSION, &major);
    cusolverGetProperty(MINOR_VERSION, &minor);
    cusolverGetProperty(PATCH_LEVEL, &patch);

    loaded = true;

    logging::info(true, "Creating CUDA handles");
    logging::up();
    logging::info(true, "Using ", gpu_properties.name);
    logging::info(true, "Memory  : ", gpu_properties.totalGlobalMem);
    logging::info(true, "CuSolver: ", major, ".", minor, ".", patch);
    logging::down();
}

void Manager::destroy_cuda() {
    if (!loaded) {
        return;
    }
    cusolverSpDestroy(handle_cusolver_sp);
    cusolverDnDestroy(handle_cusolver_dn);
    cublasDestroy_v2(handle_cublas);
    cusparseDestroy(handle_cusparse);
    loaded = false;
}

std::size_t Manager::mem_free() {
    std::size_t free_bytes = 0;
    std::size_t total_bytes = 0;
    runtime_check_cuda(cudaMemGetInfo(&free_bytes, &total_bytes));
    return free_bytes;
}

bool Manager::is_loaded() {
    return loaded;
}

Manager::~Manager() {
    destroy_cuda();
}

} // namespace cuda
} // namespace fem

#endif
