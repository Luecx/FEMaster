/**
 * @file solve_sparse_direct_cpu.cpp
 * @brief Implements CPU sparse direct solves for the FEMaster solver subsystem.
 *
 * General sparse systems are solved with LU factorization, while symmetric
 * systems use LDLT. The concrete CPU backend is selected at compile time:
 * Intel oneMKL uses PARDISO, Apple Accelerate is used for symmetric sparse
 * systems, and the portable fallback uses Eigen's native sparse solvers.
 *
 * Accelerate-enabled builds retain Eigen SparseLU for general nonsymmetric
 * matrices because Eigen only exposes Accelerate wrappers for symmetric and QR
 * factorizations. This keeps the general solve path portable across Apple SDK
 * versions while still accelerating the dominant symmetric FEM systems.
 *
 * Failed primary factorizations fall back to Eigen SparseQR so the established
 * FEMaster recovery behavior is preserved across all CPU backends.
 *
 * @author Finn Eggers
 * @date 23.08.2026
 */

#include "solve_sparse_direct.h"

#include "../../core/logging.h"
#include "../../core/timer.h"

#include <Eigen/SparseLU>
#include <Eigen/SparseQR>

#if defined(USE_MKL)
#include "../../core/config.h"
#include <Eigen/PardisoSupport>
#include <mkl.h>
#elif defined(USE_ACCELERATE)
#include <Eigen/AccelerateSupport>
#else
#include <Eigen/SparseCholesky>
#endif

namespace fem::solver::detail {

/**
 * Solves a sparse linear system on the CPU using the configured direct backend.
 *
 * General matrices are factorized with LU. Symmetric matrices are factorized
 * with LDLT so the solver can exploit their structure while still supporting
 * indefinite finite-element systems. If the selected factorization or solve
 * fails, Eigen SparseQR is used as the common recovery path.
 *
 * Accelerate-enabled builds use Eigen SparseLU for general matrices and
 * AccelerateLDLT for symmetric matrices. The latter delegates the factorization
 * and solve to Apple's sparse solver implementation while preserving Eigen's
 * sparse matrix interface.
 *
 * @param mat Sparse system matrix. Accelerate may compress its storage in place
 *            without changing its numerical entries.
 * @param rhs Dense right-hand-side matrix. Multiple columns are solved together.
 * @param matrix_type Structural matrix type used to select LU or LDLT.
 * @return Dense matrix containing the solution columns.
 */
DynamicMatrix solve_direct_cpu(SparseMatrix& mat, const DynamicMatrix& rhs, DirectSolverMatrixType matrix_type) {
    Timer t{};
    t.start();

    // General matrices require an LU factorization because no symmetry can be assumed
    if (matrix_type == DirectSolverMatrixType::General) {
        DynamicMatrix sol{};
        bool          success = false;

#if defined(USE_MKL)
        // Factorize and solve the general system with PARDISO LU
        Eigen::PardisoLU<SparseMatrix> solver{};
        solver.compute(mat);

        if (solver.info() == Eigen::Success) {
            sol     = solver.solve(rhs);
            success = solver.info() == Eigen::Success;
        }
#else
        // Keep the portable Eigen LU path for general matrices on Apple and fallback builds
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver{};
        solver.compute(mat);

        if (solver.info() == Eigen::Success) {
            sol     = solver.solve(rhs);
            success = solver.info() == Eigen::Success;
        }
#endif

        // Recover from a failed LU factorization or solve with the common QR fallback
        if (!success) {
            logging::warning(false, "General sparse LU failed; falling back to SparseQR");

            Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
            qr.compute(mat);
            sol = qr.solve(rhs);

            logging::error(qr.info() == Eigen::Success,
                "Solving general sparse system failed with SparseQR");
        }

        // Report the solve time and relative residual
        t.stop();
        logging::info(true, "Solving finished");
        logging::info(true, "Elapsed time : ", t.elapsed(), " ms");
        logging::info(true, "Residual     : ", (rhs - mat * sol).norm() / rhs.norm());
        return sol;
    }

    // Symmetric systems use LDLT so the factorization can exploit matrix symmetry
    DynamicMatrix sol{};
    bool          success = false;

#if defined(USE_MKL)
    // Configure PARDISO threading from the global FEMaster solver configuration
    mkl_set_num_threads(global_config.max_threads);
    const int mkl_max_threads = mkl_get_max_threads();

    logging::info(true, "MKL max threads          : ", mkl_max_threads);

    Eigen::PardisoLDLT<SparseMatrix> solver{};
    solver.compute(mat);

    if (solver.info() == Eigen::Success) {
        sol     = solver.solve(rhs);
        success = solver.info() == Eigen::Success;
    }
#elif defined(USE_ACCELERATE)
    // Eigen exposes Accelerate's symmetric sparse factorization directly
    mat.makeCompressed();

    Eigen::AccelerateLDLT<SparseMatrix> solver{};
    solver.compute(mat);

    if (solver.info() == Eigen::Success) {
        sol     = solver.solve(rhs);
        success = solver.info() == Eigen::Success;
    }
#else
    // Use Eigen's portable symmetric factorization when no native backend is enabled
    Eigen::SimplicialLDLT<SparseMatrix> solver{};
    solver.compute(mat);

    if (solver.info() == Eigen::Success) {
        sol     = solver.solve(rhs);
        success = solver.info() == Eigen::Success;
    }
#endif

    // Recover from a failed LDLT factorization or solve with the common QR fallback
    if (!success) {
        logging::warning(false, "Sparse LDLT failed; falling back to SparseQR");

        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr(mat);
        qr.compute(mat);
        sol = qr.solve(rhs);

        logging::error(qr.info() == Eigen::Success,
            "Solving sparse system failed with SparseQR");
    }

    // Report the solve time and relative residual
    t.stop();
    logging::info(true, "Solving finished");
    logging::info(true, "Elapsed time : ", t.elapsed(), " ms");
    logging::info(true, "Residual     : ", (rhs - mat * sol).norm() / rhs.norm());

    return sol;
}

} // namespace fem::solver::detail
