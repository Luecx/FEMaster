/**
 * @file eigen_general_buckling.cpp
 * @brief Generalized Buckling: (A - σB)^{-1} A x.
 *
 * @details Uses Spectra’s `SymGEigsShiftSolver` with:
 *            - Left operator:  y = (A - σB)^{-1} * (A * x)
 *            - Implemented via `ShiftInvertOpGeneral` (LDLT factorization)
 *
 *          A larger subspace size (ncv) is used by default: max(choose_ncv, 4k+20).
 *          Returns physical eigenvalues λ (not the transformed ν).
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 */

#include "eigen_internal.h"
#include "eigen_shift_invert_op_gen.h"

namespace fem::solver::detail {

std::vector<EigenValueVectorPair>
eigs_general_buckling(const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts)
{
    using OpA = Spectra::SparseSymMatProd<Precision>;

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = std::max(choose_ncv(n, k), std::min(n, 4 * k + 20));

    logging::info(true, "Eigen (gen): Buckling, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    ShiftInvertOpGeneral op(A, B, static_cast<Precision>(opts.sigma));
    OpA opA2(A);

    Spectra::SymGEigsShiftSolver<
        ShiftInvertOpGeneral,
        OpA,
        Spectra::GEigsMode::Buckling
    > eig(op, opA2, k, ncv_user, static_cast<Precision>(opts.sigma));

    eig.init();
    eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);

    logging::error(eig.info() == Spectra::CompInfo::Successful,
                   "Spectra generalized (Buckling) failed");

    DynamicVector evals = eig.eigenvalues(); // physical λ
    DynamicMatrix evecs = eig.eigenvectors();
    for (int i = 0; i < evals.size(); ++i)
        out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });

    t.stop();
    logging::info(true, "Eigen (gen/buckling) elapsed: ", t.elapsed(), " ms");

    return out;
}

} // namespace fem::solver::detail
