
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

    logging::info(cuda::manager.mem_free() > buffer_size,
             std::setw(60), std::left, "Allocating buffer for incomplete cholesky solve",
             std::setw(16), std::left, buffer_size,
             std::setw(16), std::left, cuda::manager.mem_free());
    logging::error(cuda::manager.mem_free() > buffer_size,
              std::setw(60), std::left, "Allocating buffer for incomplete cholesky solve",
              std::setw(16), std::left, buffer_size,
              std::setw(16), std::left, cuda::manager.mem_free());
    cuda::CudaArray<char> buffer {(size_t)buffer_size};

    logging::info(buffer_size > 1e9, "IC  Buffer storage: ", buffer_size / 1024 / 1024 / 1024.0, "Gb");

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
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#else
#endif
    logging::info(true, "");
    logging::info(true, "Solving system with N=", N, " nnz=", nnz, " using PCG (incomplete cholesky)");
    logging::up();

#ifdef SUPPORT_GPU
    if(device == GPU){
        logging::up();
        cuda::manager.create_cuda();

        // create a timer to measure time
        Timer t {};
        t.start();

        // check memory availability for precon and actual matrix on the cpu
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

        // precon
        solve_iter_gpu_precon(prec);

        SparseMatrix prec_cpu{mat};
        prec.download(prec_cpu);

        // things related to the pcg method
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
        CudaPrecision  one        = 1;
        CudaPrecision  zero       = 0;

        runtime_check_cuda(cusparseSpMV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &one, descr_A, vec_p, &zero, vec_ap, CUDA_P_TYPE,
                                                   CUSPARSE_SPMV_CSR_ALG1   , &buffer_size_ap));
        runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                                   &one, descr_L, vec_r, vec_i        , CUDA_P_TYPE,
                                                   CUSPARSE_SPSV_ALG_DEFAULT, spsv_1_descr, &buffer_size_spsv_1));
        runtime_check_cuda(cusparseSpSV_bufferSize(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    ,
                                                   &one, descr_L, vec_i, vec_z        , CUDA_P_TYPE,
                                                   CUSPARSE_SPSV_ALG_DEFAULT, spsv_2_descr, &buffer_size_spsv_2));

        logging::info(cuda::manager.mem_free() > buffer_size_ap,
                 std::setw(60), std::left, "Allocating buffer for matrix vector product",
                 std::setw(16), std::left, buffer_size_ap,
                 std::setw(16), std::left, cuda::manager.mem_free());
        logging::error(cuda::manager.mem_free() > buffer_size_ap,
                 std::setw(60), std::left, "Allocating buffer for matrix vector product",
                 std::setw(16), std::left, buffer_size_ap,
                 std::setw(16), std::left, cuda::manager.mem_free());
        cuda::CudaArray<char> buffer_ap    {buffer_size_ap};
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

        // analyse matrices for quick inversion
        runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
                                                 descr_L, vec_r, vec_i, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                                 spsv_1_descr, buffer_spsv_1));
        runtime_check_cuda(cusparseSpSV_analysis(cuda::manager.handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE    , &one,
                                                 descr_L, vec_i, vec_z, CUDA_P_TYPE, CUSPARSE_SPSV_ALG_DEFAULT,
                                                 spsv_2_descr, buffer_spsv_2));

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

        logging::info(true, "Starting iterations");
        CudaPrecision r_norm;
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

            if(k % 1000 == 0){
                logging::info("Iteration ", k, " r_norm: ", r_norm);
            }

            if(r_norm < 1e-8){
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
        logging::down();
        logging::info(true, "Running PCG method finished");
        logging::info(true, "Elapsed time: ", t.elapsed()," ms");
        logging::info(true, "iterations  : ", k);
        logging::info(true, "residual    : ", r_norm);

        DynamicVector sol{N};
        vec_x.download(sol.data());
        logging::down();
        return sol;

    } else{
#endif

        DynamicVector sol{N};
        // start the time
        Timer t {};
        t.start();

        // Convert your DenseVector into an Eigen Vector
        DynamicVector eigen_rhs = rhs;

        // Now we use Eigen's ConjugateGradient solver to solve the system
        Eigen::ConjugateGradient<Eigen::SparseMatrix<Precision>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<Precision>> cg;
        cg.compute(mat);
        cg.setTolerance(1e-12);
        logging::error(cg.info() == Eigen::Success, "Decomposition failed");
        DynamicVector eigen_sol = cg.solve(eigen_rhs);
        logging::error(cg.info() == Eigen::Success, "Solving failed");

        t.stop();
        logging::info(true, "Running PCG method finished");
        logging::info(true, "Elapsed time: ", t.elapsed()," ms");
        logging::info(true, "iterations  : ", cg.iterations());
        logging::info(true, "residual    : ", cg.error());

        // Finally, copy the result back into your rhs vector.
        for (int i = 0; i < rhs.size(); ++i) {
            sol[i] = eigen_sol(i);
        }
        logging::down();
        return sol;
#ifdef SUPPORT_GPU
    }
#endif
}
}