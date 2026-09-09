#include "cuda_vec.h"

#ifdef SUPPORT_GPU
namespace fem::cuda {
namespace {
__global__ void scalar_ratio_kernel(CudaPrecision* result,
                                    const CudaPrecision* numerator,
                                    const CudaPrecision* denominator) {
    result[0] = numerator[0] / denominator[0];
}

__global__ void scalar_ratio_pair_kernel(CudaPrecision* result,
                                         CudaPrecision* negative_result,
                                         const CudaPrecision* numerator,
                                         const CudaPrecision* denominator) {
    result[0] = numerator[0] / denominator[0];
    negative_result[0] = -result[0];
}
} // namespace

void scalar_ratio(CudaPrecision* result,
                  const CudaPrecision* numerator,
                  const CudaPrecision* denominator) {
    scalar_ratio_kernel<<<1, 1>>>(result, numerator, denominator);
    runtime_check_cuda(cudaGetLastError());
}

void scalar_ratio_pair(CudaPrecision* result,
                       CudaPrecision* negative_result,
                       const CudaPrecision* numerator,
                       const CudaPrecision* denominator) {
    scalar_ratio_pair_kernel<<<1, 1>>>(result, negative_result, numerator, denominator);
    runtime_check_cuda(cudaGetLastError());
}
} // namespace fem::cuda
#endif
