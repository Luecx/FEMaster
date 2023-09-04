
#include "solver.h"

namespace solver{

DynamicVector solve(SolverDevice device,
                    SparseMatrix& mat,
                    DynamicVector& rhs){
    runtime_assert(mat.rows() == mat.cols(), "matrix must be square");
    runtime_assert(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    runtime_assert(rhs.cols() == 1, "can only solve one equation at a time");

    const auto N   = mat.cols();
    const auto nnz = mat.nonZeros();
    DynamicVector sol{N};

#ifndef SUPPORT_GPU
    log_info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#else
#ifdef DOUBLE_PRECISION
    log_info(device != CPU, "This build does not support gpu-accelerated solving in double-precision, falling back to cpu");
    device = CPU;
#endif
#endif
    log_info(true, "");
    log_info(true, "==============================================================================");
    log_info(true, "Solving system with N=", N, " nnz=", nnz, " using cholesky decomposition");
    log_info(true, "==============================================================================");

#ifdef SUPPORT_GPU
    if (device == array::GPU) {
        // init gpu if not init yet
        cuda::manager.create_cuda();

        // Move all parts of the matrix and vector to the GPU
        mat.get_non_zero_values() >> array::GPU;
        mat.get_row_extents()     >> array::GPU;
        mat.get_col_indices()     >> array::GPU;
        sol = rhs;
        sol >> array::GPU;

        // get number of non-zero values
        int nnz = static_cast<int>(mat.get_non_zero_values().size());

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
                                                 mat.get_non_zero_values().address(array::GPU),
                                                 mat.get_row_extents().address(array::GPU),
                                                 mat.get_col_indices().address(array::GPU),
                                                 sol.address(array::GPU),
                                                 0,    // tolerance
                                                 3,    // reorder
                                                 sol.address(array::GPU),
                                                 &singularity));
        t.stop();
        log_error(singularity == -1, "decomposing system not possible");
        log_info(true, "Running PCG method finished");
        log_info(true, "Elapsed time: ", t.elapsed()," ms");

        sol >> array::CPU;

        rhs                      .free(array::GPU);
        sol                      .free(array::GPU);
        mat.get_non_zero_values().free(array::GPU);
        mat.get_row_extents()    .free(array::GPU);
        mat.get_col_indices()    .free(array::GPU);

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
        log_error(solver.info() == Eigen::Success, "Decomposition failed");
        DynamicVector eigen_sol = solver.solve(rhs);
        log_error(solver.info() == Eigen::Success, "Solving failed");

        t.stop();
        log_info(true, "Solving finished");
        log_info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        log_info(true, "residual2   : ", (rhs - mat * eigen_sol).norm() / (rhs.norm()));

        // Finally, copy the result back into your rhs vector.
        for (size_t i = 0; i < rhs.size(); ++i) {
            sol[i] = eigen_sol(i);
        }
    }
    return sol;
    //...
}

}
