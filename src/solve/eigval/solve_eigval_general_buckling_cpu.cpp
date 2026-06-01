/**
 * @file solve_eigval_general_buckling_cpu.cpp
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

#include "solve_eigval_cpu.h"
#include "solve_eigval_shiftinvert_op_general.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_general_buckling_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts)
{
    using OpA = Spectra::SparseSymMatProd<Precision>;

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = std::max(choose_ncv(n, k), std::min(n, 4 * k + 20));

    logging::info(true, "Eigval (gen): Buckling, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigvalPair> out;
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
    logging::info(true, "Eigval (gen/buckling) elapsed: ", t.elapsed(), " ms");

    return out;
}
} // namespace fem::solver::detail
