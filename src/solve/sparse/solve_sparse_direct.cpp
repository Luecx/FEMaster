#include "solve_sparse_direct.h"

#include "../../core/logging.h"
#include "../get_solver_name.h"

namespace fem::solver {

DynamicMatrix solve_direct(SolverDevice device,
                           SparseMatrix& mat,
                           const DynamicMatrix& rhs,
                           DirectSolverMatrixType matrix_type) {
    logging::error(mat.rows() == mat.cols(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    logging::error(rhs.cols() > 0, "at least one right-hand side is required");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif

    logging::info(true, "");
    logging::info(true, "Solving system with N=", mat.cols(), " nnz=", mat.nonZeros(),
                  " nrhs=", rhs.cols(),
                  " using ", get_solver_name(device, DIRECT),
                  matrix_type == DirectSolverMatrixType::SPD ? " (SPD)" : " (general)");

    logging::up();
    DynamicMatrix sol = (device == GPU)
        ? detail::solve_direct_gpu(mat, rhs, matrix_type)
        : detail::solve_direct_cpu(mat, rhs, matrix_type);
    logging::down();
    return sol;
}

} // namespace fem::solver
