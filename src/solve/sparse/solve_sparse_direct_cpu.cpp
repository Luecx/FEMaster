#include "solve_sparse_direct.h"

#include "../../core/logging.h"
#include "../../core/timer.h"
#include "../../core/config.h"

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

#ifdef USE_MKL
#include <mkl.h>
#endif

namespace fem::solver::detail {

DynamicMatrix solve_direct_cpu(SparseMatrix& mat,
                               const DynamicMatrix& rhs,
                               DirectSolverMatrixType matrix_type) {
    Timer t {};
    t.start();

    if (matrix_type == DirectSolverMatrixType::General) {
        DynamicMatrix sol;

#ifdef USE_MKL
        logging::info(true, "Using MKL PardisoLU solver");
        Eigen::PardisoLU<SparseMatrix> solver {};
        solver.compute(mat);
        if (solver.info() == Eigen::Success) {
            sol = solver.solve(rhs);
        }
#else
        logging::info(true, "Using Eigen SparseLU solver");
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver {};
        solver.compute(mat);
        if (solver.info() == Eigen::Success) {
            sol = solver.solve(rhs);
        }
#endif

        if (solver.info() != Eigen::Success) {
            logging::warning(false,
                "General sparse LU failed; falling back to SparseQR");
            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
            qr.compute(mat);
            sol = qr.solve(rhs);
            logging::error(qr.info() == Eigen::Success,
                           "Solving general sparse system failed with SparseQR");
        }

        t.stop();
        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
        logging::info(true, "residual    : ", (rhs - mat * sol).norm() / rhs.norm());
        return sol;
    }

#ifdef USE_MKL
    mkl_set_num_threads(global_config.max_threads);
    int mkl_max_threads = mkl_get_max_threads();
    logging::info(true, "MKL max threads: ", mkl_max_threads);

    logging::info(true, "Using MKL PardisoLDLT solver");
    Eigen::PardisoLDLT<SparseMatrix> solver {};

    solver.compute(mat);
    logging::warning(solver.info() == Eigen::Success, "Decomposition failed with PardisoLDLT");
    DynamicMatrix sol = solver.solve(rhs);

    if (solver.info() != Eigen::Success) {
        logging::warning(true, "Solving failed with PardisoLDLT");
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);
        logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
    }
#else
    logging::info(true, "Using Eigen SimplicialLDLT solver");

    Eigen::SimplicialLDLT<SparseMatrix> solver {};
    solver.compute(mat);
    DynamicMatrix sol;
    if (solver.info() == Eigen::Success) {
        sol = solver.solve(rhs);
    }
    if (solver.info() != Eigen::Success) {
        logging::warning(true, "SimplicialLDLT failed; falling back to SparseQR");
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);
        logging::error(qr.info() == Eigen::Success, "Solving failed with SparseQR");
    }
#endif

    t.stop();
    logging::info(true, "Solving finished");
    logging::info(true, "Elapsed time: " + std::to_string(t.elapsed()) + " ms");
    logging::info(true, "residual    : ", (rhs - mat * sol).norm() / (rhs.norm()));

    return sol;
}

} // namespace fem::solver::detail
