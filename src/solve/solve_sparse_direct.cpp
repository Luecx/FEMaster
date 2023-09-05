
#include "solver.h"

namespace solver{

DynamicVector solve_direct(SolverDevice device,
                          SparseMatrix& mat,
                          DynamicVector& rhs){
    runtime_assert(mat.rows() == mat.cols(), "matrix must be square");
    runtime_assert(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    runtime_assert(rhs.cols() == 1, "can only solve one equation at a time");

    const auto N   = mat.cols();
    const auto nnz = mat.nonZeros();
    DynamicVector sol{N};
    sol.setZero();

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#else
#ifdef CUDA_DOUBLE_PRECISION
    logging::info(device != CPU, "This build does not support gpu-accelerated solving in double-precision, falling back to cpu");
    device = CPU;
#endif
#endif
    logging::info(true, "");
    logging::info(true, "Solving system with N=", N, " nnz=", nnz, " using cholesky decomposition");

    logging::up();
#ifdef SUPPORT_GPU
    if (device == GPU) {
        // init gpu if not init yet
        cuda::manager.create_cuda();

        cuda::CudaCSR mat_gpu{mat};
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
        runtime_check_cuda(cusolverSpScsrlsvchol(cuda::manager.handle_cusolver_sp,
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
        logging::info(true, "residual2   : ", (rhs - mat * sol).norm() / (rhs.norm()));

        // destroy matrix descriptor
        runtime_check_cuda(cusparseDestroyMatDescr(descr));
    }
#endif

    //...
    if (device == CPU) {
        // start the time
        Timer t {};
        t.start();

        // Now we use Eigen's SimplicialLDLT solver to solve the system
        Eigen::SimplicialLDLT<SparseMatrix> solver {};
        solver.compute(mat);
        logging::error(solver.info() == Eigen::Success, "Decomposition failed");
        DynamicVector eigen_sol = solver.solve(rhs);
        logging::error(solver.info() == Eigen::Success, "Solving failed");

        t.stop();
        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual2   : ", (rhs - mat * eigen_sol).norm() / (rhs.norm()));

        // Finally, copy the result back into your rhs vector.
        for (size_t i = 0; i < rhs.size(); ++i) {
            sol[i] = eigen_sol(i);
        }
    }
    logging::down();
    return sol;
    //...
}

}
