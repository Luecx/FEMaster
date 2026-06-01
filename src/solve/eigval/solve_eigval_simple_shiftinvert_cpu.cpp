/**
* @file solve_eigval_simple_shiftinvert_cpu.cpp
 * @brief Simple shift–invert: A x = λ x using (A - σI)^{-1}.
 *
 * @details Uses Spectra’s `SymEigsShiftSolver` with a custom operator that
 *          applies y = (A - σI)^{-1} x. The operator uses:
 *            - MKL PardisoLDLT if compiled with USE_MKL
 *            - Eigen SimplicialLDLT otherwise
 *
 *          The operator refactorizes per call (no sigma caching).
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 */

#include "solve_eigval_cpu.h"
#include "solve_eigval_shiftinvert_op_simple.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_shiftinvert_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts)
{
    const int n   = static_cast<int>(A.rows());
    const int ncv = choose_ncv(n, k);

    logging::info(true, "Eigval (simple): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv);

    std::vector<EigvalPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    ShiftInvertOpSimple op(A, static_cast<Precision>(opts.sigma));
    Spectra::SymEigsShiftSolver<ShiftInvertOpSimple>
        eig(op, k, ncv, static_cast<Precision>(opts.sigma));

    eig.init();
    eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);

    logging::error(eig.info() == Spectra::CompInfo::Successful,
                   "Spectra simple (ShiftInvert) failed");

    DynamicVector evals = eig.eigenvalues();   // physical λ
    DynamicMatrix evecs = eig.eigenvectors();
    for (int i = 0; i < evals.size(); ++i)
        out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });

    t.stop();
    logging::info(true, "Eigval (simple/shift-invert) elapsed: ", t.elapsed(), " ms");

    return out;
}
} // namespace fem::solver::detail
