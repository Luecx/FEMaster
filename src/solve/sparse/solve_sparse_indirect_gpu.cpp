#include "solve_sparse_indirect.h"
#include "amg_preconditioner.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../cuda/assert_cuda.h"
#include "../../cuda/cuda_array.h"
#include "../../cuda/cuda_csr.h"
#include "../../cuda/cuda_defs.h"
#include "../../cuda/cuda_vec.h"

#include <algorithm>
#include <cmath>
#include <iomanip>

namespace fem::solver::detail {

DynamicMatrix solve_indirect_gpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs) {
#ifndef SUPPORT_GPU
    logging::info(true, "This build does not support gpu-accelerated solving, falling back to cpu");
    return solve_indirect_cpu(mat, rhs);
#else
    const auto N   = mat.cols();
    const auto nnz = mat.nonZeros();

    cuda::manager.create_cuda();

    Timer t {};
    t.start();

    logging::down();
    logging::info(true, "memory requirements");
    logging::up();
    logging::info(true,
             std::setw(60), std::left, "  ",
             std::setw(16), std::left, "requires",
             std::setw(16), std::left, "free");
    logging::info(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat),
             std::setw(60), std::left, "Moving sparse matrix to gpu",
             std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat),
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat),
             std::setw(59), std::left, "Moving sparse matrix to gpu",
             std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat),
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaCSR mata{mat};

    Timer amg_timer {};
    amg_timer.start();
    CpuAmgPreconditioner preconditioner;
    preconditioner.compute(mat);
    amg_timer.stop();

    logging::info(true, "CPU AMG setup finished in ", amg_timer.elapsed(), " ms");
    logging::info(true, "AMG levels: ", preconditioner.level_count());
    logging::up();
    for (std::size_t level = 0; level < preconditioner.level_count(); ++level) {
        logging::info(true,
                      "level ", level,
                      ": N=", preconditioner.level_size(level),
                      " nnz=", preconditioner.level_nnz(level));
    }
    logging::down();

    logging::info(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 5,
             std::setw(60), std::left, "Allocating vectors used for solving",
             std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 5,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 5,
             std::setw(60), std::left, "Allocating vectors used for solving",
             std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 5,
             std::setw(16), std::left, cuda::manager.mem_free());

    cuda::CudaVector vec_x {int(N)};
    cuda::CudaVector vec_r {int(N)};
    cuda::CudaVector vec_z {int(N)};
    cuda::CudaVector vec_p {int(N)};
    cuda::CudaVector vec_ap{int(N)};

    DynamicVector host_r(N);
    DynamicVector host_z(N);

    CudaPrecision val_rz;
    CudaPrecision val_pap;
    CudaPrecision val_alpha;
    CudaPrecision val_alpha2;

    cusparseSpMatDescr_t descr_A;
    runtime_check_cuda(cusparseCreateCsr(&descr_A, N, N, nnz,
                                         mata.row_ptr(),
                                         mata.col_ind(),
                                         mata.val_ptr(),
                                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_BASE_ZERO, CUDA_P_TYPE));

    size_t buffer_size_ap = 0;
    CudaPrecision one     = 1;
    CudaPrecision zero    = 0;

    runtime_check_cuda(cusparseSpMV_bufferSize(cuda::manager.handle_cusparse,
                                               CUSPARSE_OPERATION_NON_TRANSPOSE,
                                               &one, descr_A, vec_p, &zero, vec_ap,
                                               CUDA_P_TYPE, CUSPARSE_SPMV_CSR_ALG1,
                                               &buffer_size_ap));

    logging::info(cuda::manager.mem_free() > buffer_size_ap,
             std::setw(60), std::left, "Allocating buffer for matrix vector product",
             std::setw(16), std::left, buffer_size_ap,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size_ap,
             std::setw(60), std::left, "Allocating buffer for matrix vector product",
             std::setw(16), std::left, buffer_size_ap,
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer_ap{buffer_size_ap};
    logging::warning(cuda::manager.mem_free() > 1e9,
                     "Free memory is dangerously low, crashes for no reasons may occur");

    logging::info(true, "Starting iterations");
    DynamicMatrix sol = DynamicMatrix::Zero(N, rhs.cols());
    int max_iterations = 0;
    CudaPrecision max_residual = 0;

    for (Eigen::Index column = 0; column < rhs.cols(); ++column) {
        if (rhs.col(column).isZero()) {
            continue;
        }

        vec_x.clear();
        host_r = rhs.col(column);
        vec_r.upload(host_r.data());
        preconditioner.apply(host_r, host_z);
        vec_z.upload(host_z.data());
        vec_p.copy(vec_z);

        const CudaPrecision rhs_norm = static_cast<CudaPrecision>(rhs.col(column).norm());
        CudaPrecision r_norm = rhs_norm;
        int k = 0;
        for (k = 1; k <= N; ++k) {
            runtime_check_cuda(cusparseSpMV(cuda::manager.handle_cusparse,
                                            CUSPARSE_OPERATION_NON_TRANSPOSE,
                                            &one, descr_A, vec_p, &zero, vec_ap,
                                            CUDA_P_TYPE, CUSPARSE_SPMV_CSR_ALG1,
                                            buffer_ap));
            runtime_check_cuda(CUBLAS_DOT(cuda::manager.handle_cublas, N,
                                          vec_r, 1, vec_z, 1, &val_rz));
            runtime_check_cuda(CUBLAS_DOT(cuda::manager.handle_cublas, N,
                                          vec_ap, 1, vec_p, 1, &val_pap));

            logging::error(std::isfinite(val_rz) && val_rz > 0,
                           "AMG-PCG encountered non-positive r^T M^-1 r at iteration ", k,
                           ": ", val_rz);
            logging::error(std::isfinite(val_pap) && val_pap > 0,
                           "AMG-PCG encountered non-positive p^T A p at iteration ", k,
                           ": ", val_pap);

            val_alpha  = val_rz / val_pap;
            val_alpha2 = -val_alpha;

            runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, N,
                                           &val_alpha, vec_p, 1, vec_x, 1));
            runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, N,
                                           &val_alpha2, vec_ap, 1, vec_r, 1));
            runtime_check_cuda(CUBLAS_NRM(cuda::manager.handle_cublas, N,
                                          vec_r, 1, &r_norm));

            const CudaPrecision relative_residual = r_norm / rhs_norm;
            if (k <= 10 || k % 100 == 0) {
                logging::info(true,
                              "RHS ", column,
                              " iteration ", k,
                              " relative residual: ", relative_residual);
            }

            if (r_norm < 1e-8) {
                break;
            }

            vec_r.download(host_r.data());
            preconditioner.apply(host_r, host_z);
            vec_z.upload(host_z.data());

            runtime_check_cuda(CUBLAS_DOT(cuda::manager.handle_cublas, N,
                                          vec_r, 1, vec_z, 1, &val_alpha2));
            logging::error(std::isfinite(val_alpha2) && val_alpha2 > 0,
                           "AMG-PCG encountered non-positive updated r^T M^-1 r at iteration ", k,
                           ": ", val_alpha2);
            val_alpha = val_alpha2 / val_rz;

            runtime_check_cuda(CUBLAS_SCAL(cuda::manager.handle_cublas, N,
                                           &val_alpha, vec_p, 1));
            runtime_check_cuda(CUBLAS_AXPY(cuda::manager.handle_cublas, N,
                                           &one, vec_z, 1, vec_p, 1));
        }

        vec_x.download(sol.col(column).data());
        max_iterations = std::max(max_iterations, std::min(k, static_cast<int>(N)));
        max_residual = std::max(max_residual, r_norm);
    }

    runtime_check_cuda(cusparseDestroySpMat(descr_A));

    t.stop();
    logging::info(true, "Running AMG-PCG method finished");
    logging::info(true, "Elapsed time: ", t.elapsed(), " ms");
    logging::info(true, "max iterations: ", max_iterations);
    logging::info(true, "max residual  : ", max_residual);

    return sol;
#endif
}

} // namespace fem::solver::detail
