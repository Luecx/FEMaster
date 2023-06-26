
#include "solver.h"

#include <iomanip>

namespace solver{

#ifdef SUPPORT_GPU

void solve_iter_gpu_precon(cuda::CudaCSR &mat) {
    cusparseMatDescr_t mat_descr;
    csric02Info_t      info;
    int                buffer_size;

    size_t n   = mat.cols();
    size_t nnz = mat.nnz();

    runtime_check_cuda(cusparseCreateCsric02Info(&info));
    runtime_check_cuda(cusparseCreateMatDescr (&mat_descr));
    runtime_check_cuda(cusparseSetMatType     (mat_descr, CUSPARSE_MATRIX_TYPE_GENERAL));
    runtime_check_cuda(cusparseSetMatIndexBase(mat_descr, CUSPARSE_INDEX_BASE_ZERO));

    // compute buffer
    runtime_check_cuda(CUSOLV_CSRIC_BUF(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(),
                                        mat.row_ptr(), mat.col_ind(), info, &buffer_size));

    log_info(cuda::manager.mem_free() > buffer_size,
             std::setw(60), std::left, "   Allocating buffer for incomplete cholesky solve",
             std::setw(16), std::left, buffer_size,
             std::setw(16), std::left, cuda::manager.mem_free());
    log_error(cuda::manager.mem_free() > buffer_size,
              std::setw(60), std::left, "  Allocating buffer for incomplete cholesky solve",
              std::setw(16), std::left, buffer_size,
              std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer {(size_t)buffer_size};

    log_info(buffer_size > 1e9, "IC  Buffer storage: ", buffer_size / 1024 / 1024 / 1024.0, "Gb");

    // analyse and solve
    runtime_check_cuda(CUSOLV_CSRIC_ANA(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(), mat.row_ptr(),
                                        mat.col_ind(), info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));
    runtime_check_cuda(CUSOLV_CSRIC(cuda::manager.handle_cusparse, n, nnz, mat_descr, mat.val_ptr(), mat.row_ptr(),
                                    mat.col_ind(), info, CUSPARSE_SOLVE_POLICY_USE_LEVEL, buffer));

    runtime_check_cuda(cusparseDestroyCsric02Info(info));
    runtime_check_cuda(cusparseDestroyMatDescr(mat_descr));

}

#endif

DynamicVector solve_iter(SolverDevice device,
                         SparseMatrix& mat,
                         DynamicVector& rhs){
    runtime_assert(mat.cols() == mat.rows(), "matrix must be square");
    runtime_assert(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    runtime_assert(rhs.cols() == 1, "can only solve one equation at a time");

    const auto N   = mat.cols();
    const auto nnz = mat.nonZeros();

#ifndef SUPPORT_GPU
    log_info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif
    log_info(true, "");
    log_info(true, "==============================================================================");
    log_info(true, "Solving system with N=", N, " nnz=", nnz, " using PCG (incomplete cholesky)");
    log_info(true, "==============================================================================");

#ifdef SUPPORT_GPU
    if(device == GPU){
        cuda::manager.create_cuda();

        // create a timer to measure time
        Timer t {};
        t.start();

        // check memory availability for precon and actual matrix on the cpu
        log_info(true,
                 std::setw(60), std::left, "  ",
                 std::setw(16), std::left, "requires",
                 std::setw(16), std::left, "free");
        log_info(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat),
                 std::setw(60), std::left, "   Moving sparse matrix to gpu",
                 std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat),
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat),
                 std::setw(59), std::left, "  Moving sparse matrix to gpu",
                 std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat),
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaCSR prec{mat};
        log_info(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat, true),
                 std::setw(60), std::left, "   Moving preconditioned matrix to gpu",
                 std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat, true),
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > cuda::CudaCSR::estimate_mem(mat, true),
                 std::setw(59), std::left, "  Moving preconditioned matrix to gpu",
                 std::setw(16), std::left, cuda::CudaCSR::estimate_mem(mat, true),
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaCSR mata{mat, prec};

        // precon
        solve_iter_gpu_precon(prec);

        SparseMatrix prec_cpu{mat};
        prec.download(prec_cpu);

        // things related to the pcg method
        log_info(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 6,
                 std::setw(60), std::left, "   Allocating vectors used for solving",
                 std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 6,
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > cuda::CudaVector::estimate_mem(N) * 6,
                 std::setw(59), std::left, "  Allocating vectors used for solving",
                 std::setw(16), std::left, cuda::CudaVector::estimate_mem(N) * 6,
                 std::setw(16), std::left, cuda::manager.mem_free());

        cuda::CudaVector vec_x {int(N)};
        cuda::CudaVector vec_r {int(N)};
        cuda::CudaVector vec_z {int(N)};
        cuda::CudaVector vec_p {int(N)};
        cuda::CudaVector vec_ap{int(N)};
        cuda::CudaVector vec_i {int(N)};


        Precision val_rz;
        Precision val_pap;
        Precision val_alpha;
        Precision val_alpha2;

        // descriptor object for A and L
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
        cusparseSpMatSetAttribute(descr_A, CUSPARSE_SPMAT_FILL_MODE, &fill_mode, sizeof(fill_mode));
        cusparseSpMatSetAttribute(descr_L, CUSPARSE_SPMAT_FILL_MODE, &fill_mode, sizeof(fill_mode));

        // descriptors for inverse matrix operation
        cusparseSpSVDescr_t spsv_1_descr;
        cusparseSpSVDescr_t spsv_2_descr;
        cusparseSpSV_createDescr(&spsv_1_descr);
        cusparseSpSV_createDescr(&spsv_2_descr);

        // create buffers for operations
        size_t buffer_size_ap     = 0;
        size_t buffer_size_spsv_1 = 0;
        size_t buffer_size_spsv_2 = 0;
        Precision  one            = 1;
        Precision  zero           = 0;

        runtime_check_cuda(cusparseSpMV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &one, descr_A, vec_p, &zero, vec_ap, CUDA_P_TYPE,
                                                   CUSPARSE_SPMV_CSR_ALG1   , &buffer_size_ap));
        runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &one, descr_L, vec_r, vec_i        , CUDA_P_TYPE,
                                                   CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr, &buffer_size_spsv_1));
        runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    ,
                                                   &one, descr_L, vec_i, vec_z        , CUDA_P_TYPE,
                                                   CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr, &buffer_size_spsv_2));

        log_info(cuda::manager.mem_free() > buffer_size_ap,
                 std::setw(60), std::left, "   Allocating buffer for matrix vector product",
                 std::setw(16), std::left, buffer_size_ap,
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > buffer_size_ap,
                 std::setw(59), std::left, "  Allocating buffer for matrix vector product",
                 std::setw(16), std::left, buffer_size_ap,
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaArray<char> buffer_ap    {buffer_size_ap};
        log_info(cuda::manager.mem_free() > buffer_size_spsv_1,
                 std::setw(60), std::left, "   Allocating buffer 1 for triangular solve",
                 std::setw(16), std::left, buffer_size_spsv_1,
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > buffer_size_spsv_1,
                 std::setw(59), std::left, "  Allocating buffer 1 for triangular solve",
                 std::setw(16), std::left, buffer_size_spsv_1,
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaArray<char> buffer_spsv_1{buffer_size_spsv_1};
        log_info(cuda::manager.mem_free() > buffer_size_spsv_2,
                 std::setw(60), std::left, "   Allocating buffer 2 for triangular solve",
                 std::setw(16), std::left, buffer_size_spsv_2,
                 std::setw(16), std::left, cuda::manager.mem_free());
        log_error(cuda::manager.mem_free() > buffer_size_spsv_2,
                 std::setw(59), std::left, "  Allocating buffer 2 for triangular solve",
                 std::setw(16), std::left, buffer_size_spsv_2,
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaArray<char> buffer_spsv_2{buffer_size_spsv_2};
        log_warning(cuda::manager.mem_free() > 1e9, "Free memory is dangerously low, crashes for no reasons may occur");

        // analyse matrices for quick inversion
        runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                                 descr_L, vec_r, vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                                 spsv_1_descr, buffer_spsv_1));
//        buffer_spsv_1.clear();
        runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    , &one,
                                                 descr_L, vec_i, vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                                 spsv_2_descr, buffer_spsv_2));
//        buffer_spsv_2.clear();

        // beginning of iterations

        // r0 = b - A * x0  // x0 = 0 --> r0 = b
        vec_r.upload(rhs.data());
        // z0 = M^-1 * r0 --> r0 = M * z0
        runtime_check_cuda(cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                              descr_L, vec_r, vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr));
        runtime_check_cuda(cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    , &one,
                                              descr_L, vec_i, vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr));
        // p0 = z0
        vec_p.copy(vec_z);

        log_info(true, "Starting iterations");
        Precision r_norm;
        int k;
        for(k = 1; k < N; k++){
            // compute A * p
            cusparseSpMV(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, descr_A, vec_p, &zero,
                         vec_ap, CUDA_P_TYPE, CUSPARSE_SPMV_CSR_ALG1, buffer_ap);
            // compute dot product of r * z and p * Ap
            CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_r , 1, vec_z, 1, &val_rz);
            CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_ap, 1, vec_p, 1, &val_pap);

            val_alpha  = val_rz / val_pap;
            val_alpha2= -val_alpha;

            // adjust x and r
            CUBLAS_AXPY(cuda::manager.handle_cublas, N, &val_alpha , vec_p , 1, vec_x, 1);
            CUBLAS_AXPY(cuda::manager.handle_cublas, N, &val_alpha2, vec_ap, 1, vec_r, 1);

            // synchronize right here and only once.
            CUBLAS_NRM(cuda::manager.handle_cublas, N, vec_r, 1, &r_norm);

            if(r_norm < 1e-12){
                break;
            }

            // solve z = M^-1 * r
            cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one, descr_L, vec_r,
                               vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr);
            cusparseSpSV_solve(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    , &one, descr_L, vec_i,
                               vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr);

            // compute r * z again and store inside val_alpha2
            CUBLAS_DOT(cuda::manager.handle_cublas, N, vec_r, 1, vec_z, 1, &val_alpha2);

            val_alpha = val_alpha2 / val_rz;

            // update p
            // scale p = val_alpha = beta
            CUBLAS_SCAL(cuda::manager.handle_cublas, N, &val_alpha, vec_p, 1);
            CUBLAS_AXPY(cuda::manager.handle_cublas, N, &one, vec_z, 1, vec_p, 1);
        }

        // destroy things
        runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_1_descr));
        runtime_check_cuda(cusparseSpSV_destroyDescr(spsv_2_descr));
        // matrices
        runtime_check_cuda(cusparseDestroySpMat(descr_A));
        runtime_check_cuda(cusparseDestroySpMat(descr_L));

        t.stop();
        log_info(true, "Running PCG method finished");
        log_info(true, "Elapsed time: ", t.elapsed()," ms");
        log_info(true, "iterations  : ", k);
        log_info(true, "residual    : ", r_norm);

        DynamicVector sol{N};
        vec_x.download(sol.data());
        return sol;

    } else{
#endif

        DynamicVector sol{N};
//        // start the time
//        Timer t {};
//        t.start();
//
//        // First, we need to convert your ImmutableSparseMatrix into an Eigen SparseMatrix
//        Eigen::SparseMatrix<Precision> eigen_mat(N, N);
//
//        // Create a triplet list for easy construction of the sparse matrix
//        std::vector<Eigen::Triplet<Precision>> tripletList;
//        tripletList.reserve(mat.get_non_zero_values().size());
//
//        for (size_t i = 0; i < N; ++i) {
//            auto row_start = mat.get_row_extents()[i];
//            auto row_end   = mat.get_row_extents()[i + 1];
//            for (auto j = row_start; j < row_end; ++j) {
//                tripletList.emplace_back(i, mat.get_col_indices()[j], mat.get_non_zero_values()[j]);
//            }
//        }
//        eigen_mat.setFromTriplets(tripletList.begin(), tripletList.end());
//
//        // Convert your DenseVector into an Eigen Vector
//        Eigen::VectorXf eigen_rhs(N);
//        for (size_t i = 0; i < rhs.size(); ++i) {
//            eigen_rhs(i) = rhs[i];
//        }
//
//        // Now we use Eigen's SimplicialLDLT solver to solve the system
//        Eigen::ConjugateGradient<Eigen::SparseMatrix<Precision>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<Precision>> cg;
//        cg.compute(eigen_mat);
//        log_error(cg.info() == Eigen::Success, "Decomposition failed");
//        Eigen::VectorXf eigen_sol = cg.solve(eigen_rhs);
//        log_error(cg.info() == Eigen::Success, "Solving failed");
//
//        t.stop();
//        log_info(true, "Running PCG method finished");
//        log_info(true, "Elapsed time: ", t.elapsed()," ms");
//        log_info(true, "iterations  : ", cg.iterations());
//        log_info(true, "residual    : ", cg.error());
//
//        // Finally, copy the result back into your rhs vector.
//        for (size_t i = 0; i < rhs.size(); ++i) {
//            sol[i] = eigen_sol(i);
//        }
        return sol;
#ifdef SUPPORT_GPU
    }
#endif
}
}