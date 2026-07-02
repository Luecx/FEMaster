#include "solve_sparse_direct.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../cuda/assert_cuda.h"
#include "../../cuda/cuda_array.h"
#include "../../cuda/cuda_csr.h"
#include "../../cuda/cuda_defs.h"
#include "../../cuda/cuda_vec.h"

#ifdef USE_CUDSS
#include <cudss.h>
#endif

namespace fem::solver::detail {

#ifdef USE_CUDSS
namespace {
const char* cudss_status_name(cudssStatus_t status) {
    switch (status) {
        case CUDSS_STATUS_SUCCESS: return "CUDSS_STATUS_SUCCESS";
        case CUDSS_STATUS_NOT_INITIALIZED: return "CUDSS_STATUS_NOT_INITIALIZED";
        case CUDSS_STATUS_ALLOC_FAILED: return "CUDSS_STATUS_ALLOC_FAILED";
        case CUDSS_STATUS_INVALID_VALUE: return "CUDSS_STATUS_INVALID_VALUE";
        case CUDSS_STATUS_NOT_SUPPORTED: return "CUDSS_STATUS_NOT_SUPPORTED";
        case CUDSS_STATUS_EXECUTION_FAILED: return "CUDSS_STATUS_EXECUTION_FAILED";
        case CUDSS_STATUS_INTERNAL_ERROR: return "CUDSS_STATUS_INTERNAL_ERROR";
        default: return "CUDSS_STATUS_UNKNOWN";
    }
}
} // namespace
#endif

DynamicMatrix solve_direct_gpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type) {
#ifndef SUPPORT_GPU
    logging::info(true, "This build does not support gpu-accelerated solving, falling back to cpu");
    return solve_direct_cpu(mat, rhs, matrix_type);
#else
    const auto N = mat.cols();
    const auto nrhs = rhs.cols();
    DynamicMatrix sol = DynamicMatrix::Zero(N, nrhs);

    cuda::manager.create_cuda();
    cuda::CudaCSR mat_gpu{mat};
    Timer t {};

#ifdef USE_CUDSS
    logging::info(true, matrix_type == DirectSolverMatrixType::SPD
                        ? "Using cuDSS direct solver (SPD)"
                        : "Using cuDSS direct solver (general)");

    cuda::CudaArray<CudaPrecision> dense_rhs {static_cast<std::size_t>(rhs.size())};
    cuda::CudaArray<CudaPrecision> dense_sol {static_cast<std::size_t>(rhs.size())};
    dense_rhs.upload(rhs.data());

    cudssHandle_t handle {};
    cudssConfig_t config {};
    cudssData_t data {};
    cudssMatrix_t cudss_mat {};
    cudssMatrix_t cudss_rhs {};
    cudssMatrix_t cudss_sol {};

    auto runtime_check_cudss = [](cudssStatus_t status, const char* call) {
        logging::error(status == CUDSS_STATUS_SUCCESS,
                       "cuDSS call failed: ", call, " status=", cudss_status_name(status),
                       " (", static_cast<int>(status), ")");
    };

    const auto nrows = static_cast<int64_t>(N);
    const auto ncols = static_cast<int64_t>(N);
    const auto nnz_gpu = static_cast<int64_t>(mat_gpu.nnz());

    int* row_offsets = mat_gpu.row_ptr();
    int* row_start = row_offsets;
    int* col_indices = mat_gpu.col_ind();
    CudaPrecision* values = mat_gpu.val_ptr();
    CudaPrecision* rhs_values = dense_rhs;
    CudaPrecision* sol_values = dense_sol;
    const cudssMatrixType_t cudss_matrix_type =
        (matrix_type == DirectSolverMatrixType::SPD) ? CUDSS_MTYPE_SPD : CUDSS_MTYPE_GENERAL;
    const cudssMatrixViewType_t cudss_matrix_view =
        (matrix_type == DirectSolverMatrixType::SPD) ? CUDSS_MVIEW_UPPER : CUDSS_MVIEW_FULL;

    t.start();

    runtime_check_cudss(cudssCreate(&handle), "cudssCreate");
    runtime_check_cudss(cudssConfigCreate(&config), "cudssConfigCreate");
    runtime_check_cudss(cudssDataCreate(handle, &data), "cudssDataCreate");
    runtime_check_cudss(cudssMatrixCreateCsr(&cudss_mat,
                                             nrows,
                                             ncols,
                                             nnz_gpu,
                                             row_start,
                                             nullptr,
                                             col_indices,
                                             values,
                                             CUDA_R_32I,
                                             CUDA_P_TYPE,
                                             cudss_matrix_type,
                                             cudss_matrix_view,
                                             CUDSS_BASE_ZERO),
                        "cudssMatrixCreateCsr");
    runtime_check_cudss(cudssMatrixCreateDn(&cudss_rhs,
                                            nrows,
                                            static_cast<int64_t>(nrhs),
                                            nrows,
                                            rhs_values,
                                            CUDA_P_TYPE,
                                            CUDSS_LAYOUT_COL_MAJOR),
                        "cudssMatrixCreateDn(rhs)");
    runtime_check_cudss(cudssMatrixCreateDn(&cudss_sol,
                                            nrows,
                                            static_cast<int64_t>(nrhs),
                                            nrows,
                                            sol_values,
                                            CUDA_P_TYPE,
                                            CUDSS_LAYOUT_COL_MAJOR),
                        "cudssMatrixCreateDn(sol)");

    runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, config, data, cudss_mat, cudss_sol, cudss_rhs),
                        "cudssExecute(analysis)");
    runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, config, data, cudss_mat, cudss_sol, cudss_rhs),
                        "cudssExecute(factorization)");
    runtime_check_cudss(cudssExecute(handle, CUDSS_PHASE_SOLVE, config, data, cudss_mat, cudss_sol, cudss_rhs),
                        "cudssExecute(solve)");

    t.stop();
    dense_sol.download(sol.data());

    runtime_check_cudss(cudssMatrixDestroy(cudss_mat), "cudssMatrixDestroy(mat)");
    runtime_check_cudss(cudssMatrixDestroy(cudss_sol), "cudssMatrixDestroy(sol)");
    runtime_check_cudss(cudssMatrixDestroy(cudss_rhs), "cudssMatrixDestroy(rhs)");
    runtime_check_cudss(cudssDataDestroy(handle, data), "cudssDataDestroy");
    runtime_check_cudss(cudssConfigDestroy(config), "cudssConfigDestroy");
    runtime_check_cudss(cudssDestroy(handle), "cudssDestroy");
#else
    logging::error(matrix_type == DirectSolverMatrixType::SPD,
                   "GPU direct solve for non-SPD systems requires cuDSS; "
                   "cuSolver sparse direct path only supports Cholesky/SPD here");
    logging::info(nrhs > 1,
                  "Legacy cuSolver solves multiple right-hand sides sequentially; "
                  "enable cuDSS for a native multi-RHS solve");

    cuda::CudaVector vec_rhs {int(N)};
    cuda::CudaVector vec_sol {int(N)};

    int nnz = static_cast<int>(mat_gpu.nnz());

    cusparseMatDescr_t descr;
    runtime_check_cuda(cusparseCreateMatDescr(&descr));
    runtime_check_cuda(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    runtime_check_cuda(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

    t.start();
    for (Eigen::Index column = 0; column < nrhs; ++column) {
        vec_rhs.upload(rhs.col(column).data());
        int singularity;
        runtime_check_cuda(CUSOLV_CHOLESKY(cuda::manager.handle_cusolver_sp,
                                           N,
                                           nnz,
                                           descr,
                                           mat_gpu.val_ptr(),
                                           mat_gpu.row_ptr(),
                                           mat_gpu.col_ind(),
                                           vec_rhs,
                                           0,
                                           3,
                                           vec_sol,
                                           &singularity));
        logging::error(singularity == -1, "decomposing system not possible for RHS ", column);
        vec_sol.download(sol.col(column).data());
    }
    t.stop();

    runtime_check_cuda(cusparseDestroyMatDescr(descr));
#endif

    logging::info(true, "Solving finished");
    logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    logging::info(true, "residual    : ", (rhs - mat * sol).norm() / (rhs.norm()));
    return sol;
#endif
}

} // namespace fem::solver::detail
