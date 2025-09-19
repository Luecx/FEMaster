/******************************************************************************
 * @file eigen_general_shiftinvert.cpp
 * @brief Generalized ShiftInvert: (A - σB)^{-1} B x.
 *
 * @details Uses Spectra’s `SymGEigsShiftSolver` with:
 *            - Left operator:  y = (A - σB)^{-1} * (B * x)
 *            - Implemented via `ShiftInvertOpGeneral` (LDLT factorization)
 *
 *          Subspace size is chosen automatically. Returns physical eigenvalues.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 ******************************************************************************/

#include "eigen_internal.h"
#include "eigen_shift_invert_op_gen.h"

namespace fem::solver::detail {

std::vector<EigenValueVectorPair>
eigs_general_shiftinvert(const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts)
{
    using OpB = Spectra::SparseSymMatProd<Precision>;

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = choose_ncv(n, k);

    logging::info(true, "Eigen (gen): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigenValueVectorPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    ShiftInvertOpGeneral op(A, B, static_cast<Precision>(opts.sigma));
    OpB opB(B);

    Spectra::SymGEigsShiftSolver<
        ShiftInvertOpGeneral,
        OpB,
        Spectra::GEigsMode::ShiftInvert
    > eig(op, opB, k, ncv_user, static_cast<Precision>(opts.sigma));

    eig.init();
    eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);

    logging::error(eig.info() == Spectra::CompInfo::Successful,
                   "Spectra generalized (ShiftInvert) failed");

    DynamicVector evals = eig.eigenvalues();
    DynamicMatrix evecs = eig.eigenvectors();
    for (int i = 0; i < evals.size(); ++i)
        out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });

    t.stop();
    logging::info(true, "Eigen (gen/shift-invert) elapsed: ", t.elapsed(), " ms");

    return out;
}

} // namespace fem::solver::detail
