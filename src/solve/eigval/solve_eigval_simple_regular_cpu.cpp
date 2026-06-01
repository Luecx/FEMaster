/**
* @file solve_eigval_simple_regular_cpu.cpp
 * @brief Simple regular eigval solver: A x = λ x (A symmetric).
 *
 * @details Uses Spectra’s `SymEigsSolver` with the matrix–vector operator
 *          `SparseSymMatProd`. No inner linear solves are performed.
 *
 *          The subspace size `ncv` is chosen automatically from k and n.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 */

#include "solve_eigval_cpu.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_regular_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts)
{
    using OpA = Spectra::SparseSymMatProd<Precision>;

    const int n   = static_cast<int>(A.rows());
    const int ncv = choose_ncv(n, k);

    logging::info(true, "Eigval (simple): Regular, ncv=", ncv);

    std::vector<EigvalPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    OpA opA(A);
    Spectra::SymEigsSolver<OpA> eig(opA, k, ncv);

    eig.init();
    eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);

    logging::error(eig.info() == Spectra::CompInfo::Successful,
                   "Spectra simple (Regular) failed");

    DynamicVector evals = eig.eigenvalues();
    DynamicMatrix evecs = eig.eigenvectors();
    for (int i = 0; i < evals.size(); ++i)
        out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });

    t.stop();
    logging::info(true, "Eigval (simple/regular) elapsed: ", t.elapsed(), " ms");

    return out;
}
} // namespace fem::solver::detail
