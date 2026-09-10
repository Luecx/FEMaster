/**
 * @file cuda_defs.h
 * @brief Maps precision-specific CUDA routines to unified identifiers.
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#ifdef CUDA_DOUBLE_PRECISION
#define CUDA_P_TYPE      CUDA_R_64F
#define CUBLAS_DOT       cublasDdot
#define CUBLAS_NRM       cublasDnrm2
#define CUBLAS_AXPY      cublasDaxpy
#define CUBLAS_SCAL      cublasDscal
#define CUSOLV_CSRIC_BUF cusparseDcsric02_bufferSize
#define CUSOLV_CSRIC_ANA cusparseDcsric02_analysis
#define CUSOLV_CSRIC     cusparseDcsric02
#define CUSOLV_CHOLESKY  cusolverSpDcsrlsvchol
#else
#define CUDA_P_TYPE      CUDA_R_32F
#define CUBLAS_DOT       cublasSdot
#define CUBLAS_NRM       cublasSnrm2
#define CUBLAS_AXPY      cublasSaxpy
#define CUBLAS_SCAL      cublasSscal
#define CUSOLV_CSRIC_BUF cusparseScsric02_bufferSize
#define CUSOLV_CSRIC_ANA cusparseScsric02_analysis
#define CUSOLV_CSRIC     cusparseScsric02
#define CUSOLV_CHOLESKY  cusolverSpScsrlsvchol
#endif

#ifdef USE_CUDSS
#include <cudss.h>

namespace fem::cuda::detail {

#ifdef CUDA_DOUBLE_PRECISION
inline constexpr cudssDataType_t cudss_precision_type = CUDSS_R_64F;
#else
inline constexpr cudssDataType_t cudss_precision_type = CUDSS_R_32F;
#endif

// cuDSS 0.8 uses its own data-type enum and separate offset/index types.
template<typename LegacyIndexType, typename LegacyValueType>
inline cudssStatus_t cudss_matrix_create_csr_compat(cudssMatrix_t* matrix,
                                                     int64_t nrows,
                                                     int64_t ncols,
                                                     int64_t nnz,
                                                     const void* row_start,
                                                     const void* row_end,
                                                     const void* col_indices,
                                                     const void* values,
                                                     LegacyIndexType,
                                                     LegacyValueType,
                                                     cudssMatrixType_t matrix_type,
                                                     cudssMatrixViewType_t matrix_view,
                                                     cudssIndexBase_t index_base) {
    return ::cudssMatrixCreateCsr(matrix,
                                  nrows,
                                  ncols,
                                  nnz,
                                  row_start,
                                  row_end,
                                  col_indices,
                                  values,
                                  CUDSS_R_32I,
                                  CUDSS_R_32I,
                                  cudss_precision_type,
                                  matrix_type,
                                  matrix_view,
                                  index_base);
}

template<typename LegacyValueType>
inline cudssStatus_t cudss_matrix_create_dn_compat(cudssMatrix_t* matrix,
                                                    int64_t nrows,
                                                    int64_t ncols,
                                                    int64_t ld,
                                                    const void* values,
                                                    LegacyValueType,
                                                    cudssLayout_t layout) {
    return ::cudssMatrixCreateDn(matrix,
                                 nrows,
                                 ncols,
                                 ld,
                                 values,
                                 cudss_precision_type,
                                 layout);
}

} // namespace fem::cuda::detail

#define cudssMatrixCreateCsr(...) ::fem::cuda::detail::cudss_matrix_create_csr_compat(__VA_ARGS__)
#define cudssMatrixCreateDn(...)  ::fem::cuda::detail::cudss_matrix_create_dn_compat(__VA_ARGS__)
#endif
