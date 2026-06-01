#include "solve_sparse_direct.h"

#include "../../core/logging.h"

namespace fem::solver {

DynamicVector solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           DynamicVector& rhs,
                           DirectSolverMatrixType matrix_type) {
    logging::error(mat.rows() == mat.cols(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    logging::error(rhs.cols() == 1, "can only solve one equation at a time");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif

    logging::info(true, "");
    logging::info(true, "Solving system with N=", mat.cols(), " nnz=", mat.nonZeros(),
                  matrix_type == DirectSolverMatrixType::SPD
                      ? " using cholesky decomposition"
                      : " using general direct decomposition");

    logging::up();
    DynamicVector sol = (device == GPU)
        ? detail::solve_direct_gpu(mat, rhs, matrix_type)
        : detail::solve_direct_cpu(mat, rhs, matrix_type);
    logging::down();
    return sol;
}

} // namespace fem::solver
