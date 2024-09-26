
#include "solver.h"

#ifdef USE_MKL
#include <mkl.h>
#endif

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
#endif
    logging::info(true, "");
    logging::info(true, "Solving system with N=", N, " nnz=", nnz, " using cholesky decomposition");

    logging::up();
#ifdef SUPPORT_GPU
    if (device == GPU) {
        // init gpu if not init yet
        std::cout << "trying to create cuda" << std::endl;
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
        logging::error(solver.info() == Eigen::Success, "Decomposition failed with PardisoLDLT");
        DynamicVector eigen_sol = solver.solve(rhs);
        logging::error(solver.info() == Eigen::Success, "Solving failed with PardisoLDLT");

#else
        // Use Eigen's SimplicialLDLT solver
        logging::info(true, "Using Eigen SimplicialLDLT solver");

        Eigen::SimplicialLDLT<SparseMatrix> solver {};
        solver.compute(mat);
        logging::error(solver.info() == Eigen::Success, "Decomposition failed with SimplicialLDLT");
        DynamicVector eigen_sol = solver.solve(rhs);
        logging::error(solver.info() == Eigen::Success, "Solving failed with SimplicialLDLT");

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
    //...
}

}
