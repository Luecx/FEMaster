#include "solve_sparse_indirect.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/IterativeLinearSolvers>

#include <algorithm>

namespace fem::solver::detail {

DynamicMatrix solve_indirect_cpu(SparseMatrix& mat,
                                 const DynamicMatrix& rhs) {
    Timer t {};
    t.start();

    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<Precision>,
        Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<Precision>
    > cg;
    cg.compute(mat);
    cg.setTolerance(1e-12);
    logging::error(cg.info() == Eigen::Success, "Decomposition failed");

    DynamicMatrix sol(mat.rows(), rhs.cols());
    Eigen::Index max_iterations = 0;
    Precision max_error = 0;
    for (Eigen::Index column = 0; column < rhs.cols(); ++column) {
        sol.col(column) = cg.solve(rhs.col(column));
        logging::error(cg.info() == Eigen::Success, "Solving failed for RHS ", column);
        max_iterations = std::max(max_iterations, cg.iterations());
        max_error = std::max(max_error, cg.error());
    }

    t.stop();
    logging::info(true, "Running PCG method finished");
    logging::info(true, "Elapsed time: ", t.elapsed(), " ms");
    logging::info(true, "max iterations: ", max_iterations);
    logging::info(true, "max residual  : ", max_error);

    return sol;
}

} // namespace fem::solver::detail
