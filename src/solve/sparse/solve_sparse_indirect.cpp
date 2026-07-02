#include "solve_sparse_indirect.h"

#include "../../core/logging.h"
#include "../get_solver_name.h"

namespace fem::solver {

DynamicMatrix solve_indirect(SolverDevice device,
                             SparseMatrix& mat,
                             const DynamicMatrix& rhs) {
    logging::error(mat.cols() == mat.rows(), "matrix must be square");
    logging::error(rhs.rows() == mat.rows(), "missmatch of rhs and matrix");
    logging::error(rhs.cols() > 0, "at least one right-hand side is required");

#ifndef SUPPORT_GPU
    logging::info(device != CPU, "This build does not support gpu-accelerated solving, falling back to cpu");
    device = CPU;
#endif

    logging::info(true, "");
    logging::info(true, "Solving system with N=", mat.cols(), " nnz=", mat.nonZeros(),
                  " nrhs=", rhs.cols(),
                  " using ", get_solver_name(device, INDIRECT));
    logging::up();

    DynamicMatrix sol = (device == GPU)
        ? detail::solve_indirect_gpu(mat, rhs)
        : detail::solve_indirect_cpu(mat, rhs);

    logging::down();
    return sol;
}

} // namespace fem::solver
