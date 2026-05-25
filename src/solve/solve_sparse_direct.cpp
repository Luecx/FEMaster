
#include "solver.h"

#include "../core/logging.h"
#include "../core/timer.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>

#ifdef USE_MKL
#include <mkl.h>
#endif

#ifdef USE_CUDSS
#include <cudss.h>
#endif

namespace fem::solver{
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
}
#endif

DynamicVector solve_direct(SolverDevice device,
                          SparseMatrix& mat,
                          DynamicVector& rhs) {
    logging::error(mat.rows() == mat.cols(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    logging::error(rhs.cols() == 1, "can only solve one equation at a time");

    const auto N   = mat.cols();
    const auto nnz = mat.nonZeros();
    DynamicVector sol{N};
    sol.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#else
#endif
    logging::info(true, "");
    logging::info(true, "Solving system with N=", N, " nnz=", nnz, " using cholesky decomposition");

    logging::up();
#ifdef SUPPORT_GPU
    if (device == GPU) {
        // init gpu if not init yet
        cuda::manager.create_cuda();

        cuda::CudaCSR mat_gpu{mat};

#ifdef USE_CUDSS
        logging::info(true, "Using cuDSS direct solver");

        cuda::CudaArray<CudaPrecision> vec_rhs {static_cast<std::size_t>(rhs.size())};
        cuda::CudaArray<CudaPrecision> vec_sol {static_cast<std::size_t>(rhs.size())};
        vec_rhs.upload(rhs.data());

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
        CudaPrecision* rhs_values = vec_rhs;
        CudaPrecision* sol_values = vec_sol;

        Timer t {};
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
                                                 CUDSS_MTYPE_SPD,
                                                 CUDSS_MVIEW_UPPER,
                                                 CUDSS_BASE_ZERO),
                            "cudssMatrixCreateCsr");
        runtime_check_cudss(cudssMatrixCreateDn(&cudss_rhs,
                                                nrows,
                                                1,
                                                nrows,
                                                rhs_values,
                                                CUDA_P_TYPE,
                                                CUDSS_LAYOUT_COL_MAJOR),
                            "cudssMatrixCreateDn(rhs)");
        runtime_check_cudss(cudssMatrixCreateDn(&cudss_sol,
                                                nrows,
                                                1,
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
        vec_sol.download(sol.data());

        runtime_check_cudss(cudssMatrixDestroy(cudss_mat), "cudssMatrixDestroy(mat)");
        runtime_check_cudss(cudssMatrixDestroy(cudss_sol), "cudssMatrixDestroy(sol)");
        runtime_check_cudss(cudssMatrixDestroy(cudss_rhs), "cudssMatrixDestroy(rhs)");
        runtime_check_cudss(cudssDataDestroy(handle, data), "cudssDataDestroy");
        runtime_check_cudss(cudssConfigDestroy(config), "cudssConfigDestroy");
        runtime_check_cudss(cudssDestroy(handle), "cudssDestroy");

        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual    : ", (rhs - mat * sol).norm() / (rhs.norm()));
#else
        cuda::CudaVector vec_rhs {int(rhs.size())};
        cuda::CudaVector vec_sol {int(rhs.size())};
        vec_rhs.upload(rhs.data());
        vec_sol.upload(rhs.data());

        // get number of non-zero values
        int nnz = static_cast<int>(mat_gpu.nnz());

        // create matrix descriptor
        cusparseMatDescr_t descr;
        runtime_check_cuda(cusparseCreateMatDescr(&descr));
        runtime_check_cuda(cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL));
        runtime_check_cuda(cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO));

        // stop the time
        Timer t {};
        t.start();
        int singularity;
        runtime_check_cuda(CUSOLV_CHOLESKY(cuda::manager.handle_cusolver_sp,
                                                 N,
                                                 nnz,
                                                 descr,
                                                 mat_gpu.val_ptr(),
                                                 mat_gpu.row_ptr(),
                                                 mat_gpu.col_ind(),
                                                 vec_rhs,
                                                 0,    // tolerance
                                                 3,    // reorder
                                                 vec_sol,
                                                 &singularity));
        t.stop();
        logging::error(singularity == -1, "decomposing system not possible");
        vec_sol.download(sol.data());

        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual    : ", (rhs - mat * sol).norm() / (rhs.norm()));

        // destroy matrix descriptor
        runtime_check_cuda(cusparseDestroyMatDescr(descr));
#endif
    }
#endif

    //...
    if (device == CPU) {
        // Start the time
        Timer t {};
        t.start();

#ifdef USE_MKL
        mkl_set_num_threads(global_config.max_threads);
        int mkl_max_threads = mkl_get_max_threads();
        logging::info(true, "MKL max threads: ", mkl_max_threads);

        // Use MKL PardisoLDLT solver
        logging::info(true, "Using MKL PardisoLDLT solver");
        Eigen::PardisoLDLT<SparseMatrix> solver {};

        solver.compute(mat);
        logging::warning(solver.info() == Eigen::Success, "Decomposition failed with PardisoLDLT");
        DynamicVector eigen_sol = solver.solve(rhs);

        if (solver.info() != Eigen::Success) {
            logging::warning(true, "Solving failed with PardisoLDLT");
            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
            qr.compute(mat);
            eigen_sol = qr.solve(rhs);
            logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
        }
#else
        // Use Eigen's SimplicialLDLT solver
        logging::info(true, "Using Eigen SimplicialLDLT solver");

        Eigen::SimplicialLDLT<SparseMatrix> solver {};
        solver.compute(mat);
        DynamicVector eigen_sol;
        if (solver.info() == Eigen::Success) {
            eigen_sol = solver.solve(rhs);
        }
        if (solver.info() != Eigen::Success) {
            logging::warning(true, "SimplicialLDLT failed; falling back to SparseQR");
            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
            qr.compute(mat);
            eigen_sol = qr.solve(rhs);
            logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
        }

#endif

        t.stop();
        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual    : ", (rhs - mat * eigen_sol).norm() / (rhs.norm()));

        // Finally, copy the result back into your sol vector.
        for (int i = 0; i < rhs.size(); ++i) {
            sol[i] = eigen_sol(i);
        }
    }

    logging::down();
    return sol;
}
}
