#include "solve_sparse_indirect.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../cuda/assert_cuda.h"
#include "../../cuda/cuda_array.h"
#include "../../cuda/cuda_csr.h"
#include "../../cuda/cuda_defs.h"
#include "../../cuda/cuda_vec.h"

#include <iomanip>

namespace fem::solver::detail {

#ifdef SUPPORT_GPU
namespace {
void solve_indirect_gpu_precon(cuda::CudaCSR& mat) {
    cusparseMatDescr_t mat_descr;
    csric02Info_t      info;
    int                buffer_size;

    size_t n   = mat.cols();
    size_t nnz = mat.nnz();

    runtime_check_cuda(cusparseCreateCsric02Info(&info));
    runtime_check_cuda(cusparseCreateMatDescr (&mat_descr));
    runtime_check_cuda(cusparseSetMatType     (mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    runtime_check_cuda(cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO));

    runtime_check_cuda(CUSOLV_CSRIC_BUF(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(),
                                        mat.row_ptr(), mat.col_ind(), info, &buffer_size));

    logging::info(cuda::manager.mem_free() > buffer_size,
             std::setw(60), std::left, "Allocating buffer for incomplete cholesky solve",
             std::setw(16), std::left, buffer_size,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size,
              std::setw(60), std::left, "Allocating buffer for incomplete cholesky solve",
              std::setw(16), std::left, buffer_size,
              std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer {static_cast<size_t>(buffer_size)};

    logging::info(buffer_size > 1e9, "IC  Buffer storage: ", buffer_size / 1024 / 1024 / 1024.0, "Gb");

    runtime_check_cuda(CUSOLV_CSRIC_ANA(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(), mat.row_ptr(),
                                        mat.col_ind(), info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
    runtime_check_cuda(CUSOLV_CSRIC(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(), mat.row_ptr(),
                                    mat.col_ind(), info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));

    runtime_check_cuda(cusparseDestroyCsric02Info(info));
    runtime_check_cuda(cusparseDestroyMatDescr(mat_descr));
}
} // namespace
#endif

DynamicVector solve_indirect_gpu(SparseMatrix& mat,
                                 DynamicVector& rhs) {
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
    cuda::CudaCSR prec{mat};
    logging::info(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat, true),
             std::setw(60), std::left, "Moving preconditioned matrix to gpu",
             std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat, true),
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat, true),
             std::setw(60), std::left, "Moving preconditioned matrix to gpu",
             std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat, true),
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaCSR mata{mat, prec};

    solve_indirect_gpu_precon(prec);

    SparseMatrix prec_cpu{mat};
    prec.download(prec_cpu);

    logging::info(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 6,
             std::setw(60), std::left, "Allocating vectors used for solving",
             std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 6,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 6,
             std::setw(60), std::left, "Allocating vectors used for solving",
             std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 6,
             std::setw(16), std::left, cuda::manager.mem_free());

    cuda::CudaVector vec_x {int(N)};
    cuda::CudaVector vec_r {int(N)};
    cuda::CudaVector vec_z {int(N)};
    cuda::CudaVector vec_p {int(N)};
    cuda::CudaVector vec_ap{int(N)};
    cuda::CudaVector vec_i {int(N)};

    CudaPrecision val_rz;
    CudaPrecision val_pap;
    CudaPrecision val_alpha;
    CudaPrecision val_alpha2;

    cusparseSpMatDescr_t descr_A;
    cusparseSpMatDescr_t descr_L;
    runtime_check_cuda(cusparseCreateCsr(&descr_A, N, N, nnz,
                                         mata.row_ptr(),
                                         mata.col_ind(),
                                         mata.val_ptr(),
                                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_BASE_ZERO, CUDA_P_TYPE));
    runtime_check_cuda(cusparseCreateCsr(&descr_L, N, N, nnz,
                                         prec.row_ptr(),
                                         prec.col_ind(),
                                         prec.val_ptr(),
                                         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                         CUSPARSE_INDEX_BASE_ZERO, CUDA_P_TYPE));
    auto fill_mode = CUSPARSE_FILL_MODE_LOWER;
    runtime_check_cuda(cusparseSpMatSetAttribute(descr_L, CUSPARSE_SPMAT_FILL_MODE, &fill_mode, sizeof(fill_mode)));

    cusparseSpSVDescr_t spsv_1_descr;
    cusparseSpSVDescr_t spsv_2_descr;
    cusparseSpSV_createDescr(&spsv_1_descr);
    cusparseSpSV_createDescr(&spsv_2_descr);

    size_t buffer_size_ap     = 0;
    size_t buffer_size_spsv_1 = 0;
    size_t buffer_size_spsv_2 = 0;
    CudaPrecision one         = 1;
    CudaPrecision zero        = 0;

    runtime_check_cuda(cusparseSpMV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                               &one, descr_A, vec_p, &zero, vec_ap, CUDA_P_TYPE,
                                               CUSPARSE_SPMV_CSR_ALG1, &buffer_size_ap));
    runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                               &one, descr_L, vec_r, vec_i, CUDA_P_TYPE,
                                               CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr, &buffer_size_spsv_1));
    runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE,
                                               &one, descr_L, vec_i, vec_z, CUDA_P_TYPE,
                                               CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr, &buffer_size_spsv_2));

    logging::info(cuda::manager.mem_free() > buffer_size_ap,
             std::setw(60), std::left, "Allocating buffer for matrix vector product",
             std::setw(16), std::left, buffer_size_ap,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size_ap,
             std::setw(60), std::left, "Allocating buffer for matrix vector product",
             std::setw(16), std::left, buffer_size_ap,
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer_ap{buffer_size_ap};
    logging::info(cuda::manager.mem_free() > buffer_size_spsv_1,
             std::setw(60), std::left, "Allocating buffer 1 for triangular solve",
             std::setw(16), std::left, buffer_size_spsv_1,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size_spsv_1,
             std::setw(60), std::left, "Allocating buffer 1 for triangular solve",
             std::setw(16), std::left, buffer_size_spsv_1,
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer_spsv_1{buffer_size_spsv_1};
    logging::info(cuda::manager.mem_free() > buffer_size_spsv_2,
             std::setw(60), std::left, "Allocating buffer 2 for triangular solve",
             std::setw(16), std::left, buffer_size_spsv_2,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size_spsv_2,
             std::setw(60), std::left, "Allocating buffer 2 for triangular solve",
             std::setw(16), std::left, buffer_size_spsv_2,
             std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer_spsv_2{buffer_size_spsv_2};
    logging::warning(cuda::manager.mem_free() > 1e9, "Free memory is dangerously low, crashes for no reasons may occur");

    runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                             descr_L, vec_r, vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                             spsv_1_descr, buffer_spsv_1));
    runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE, &one,
                                             descr_L, vec_i, vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                             spsv_2_descr, buffer_spsv_2));

    vec_r.upload(rhs.data());
    runtime_check_cuda(cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                          descr_L, vec_r, vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr));
    runtime_check_cuda(cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE, &one,
                                          descr_L, vec_i, vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr));
    vec_p.copy(vec_z);

    logging::info(true, "Starting iterations");
    CudaPrecision r_norm;
    int k;
    for (k = 1; k < N; k++) {
        cusparseSpMV(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, descr_A, vec_p, &zero,
                     vec_ap, CUDA_P_TYPE, CUSPARSE_SPMV_CSR_ALG1, buffer_ap);
        CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_r, 1, vec_z, 1, &val_rz);
        CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_ap, 1, vec_p, 1, &val_pap);

        val_alpha  = val_rz / val_pap;
        val_alpha2 = -val_alpha;

        CUBLAS_AXPY(cuda::manager.handle_cublas, N, &val_alpha, vec_p, 1, vec_x, 1);
        CUBLAS_AXPY(cuda::manager.handle_cublas, N, &val_alpha2, vec_ap, 1, vec_r, 1);

        CUBLAS_NRM(cuda::manager.handle_cublas, N, vec_r, 1, &r_norm);

        if (k % 1000 == 0) {
            logging::info("Iteration ", k, " r_norm: ", r_norm);
        }

        if (r_norm < 1e-8) {
            break;
        }

        cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, descr_L, vec_r,
                           vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr);
        cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE, &one, descr_L, vec_i,
                           vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr);

        CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_r, 1, vec_z, 1, &val_alpha2);

        val_alpha = val_alpha2 / val_rz;

        CUBLAS_SCAL(cuda::manager.handle_cublas, N, &val_alpha, vec_p, 1);
        CUBLAS_AXPY(cuda::manager.handle_cublas, N, &one, vec_z, 1, vec_p, 1);
    }

    runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_1_descr));
    runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_2_descr));
    runtime_check_cuda(cusparseDestroySpMat(descr_A));
    runtime_check_cuda(cusparseDestroySpMat(descr_L));

    t.stop();
    logging::info(true, "Running PCG method finished");
    logging::info(true, "Elapsed time: ", t.elapsed(), " ms");
    logging::info(true, "iterations  : ", k);
    logging::info(true, "residual    : ", r_norm);

    DynamicVector sol{N};
    vec_x.download(sol.data());
    return sol;
#endif
}

} // namespace fem::solver::detail
