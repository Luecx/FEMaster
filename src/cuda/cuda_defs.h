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
