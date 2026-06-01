#include "solve_sparse_indirect.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/IterativeLinearSolvers>

namespace fem::solver::detail {

DynamicVector solve_indirect_cpu(SparseMatrix& mat,
                                 DynamicVector& rhs) {
    Timer t {};
    t.start();

    DynamicVector indirect_rhs = rhs;

    Eigen::ConjugateGradient<
        Eigen::SparseMatrix<Precision>,
        Eigen::Lower | Eigen::Upper,
        Eigen::IncompleteCholesky<Precision>
    > cg;
    cg.compute(mat);
    cg.setTolerance(1e-12);
    logging::error(cg.info() == Eigen::Success, "Decomposition failed");
    DynamicVector sol = cg.solve(indirect_rhs);
    logging::error(cg.info() == Eigen::Success, "Solving failed");

    t.stop();
    logging::info(true, "Running PCG method finished");
    logging::info(true, "Elapsed time: ", t.elapsed(), " ms");
    logging::info(true, "iterations  : ", cg.iterations());
    logging::info(true, "residual    : ", cg.error());

    return sol;
}

} // namespace fem::solver::detail
