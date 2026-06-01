#include "solve_sparse_indirect.h"

#include "../../core/logging.h"

namespace fem::solver {

DynamicVector solve_indirect(SolverDevice device,
                             SparseMatrix& mat,
                             DynamicVector& rhs) {
    logging::error(mat.cols() == mat.rows(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    logging::error(rhs.cols() == 1, "can only solve one equation at a time");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif

    logging::info(true, "");
    logging::info(true, "Solving system with N=", mat.cols(), " nnz=", mat.nonZeros(),
                  " using PCG (incomplete cholesky)");
    logging::up();

    DynamicVector sol = (device == GPU)
        ? detail::solve_indirect_gpu(mat, rhs)
        : detail::solve_indirect_cpu(mat, rhs);

    logging::down();
    return sol;
}

} // namespace fem::solver
