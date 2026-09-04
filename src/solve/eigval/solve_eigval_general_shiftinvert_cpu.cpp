/**
 * @file solve_eigval_general_shiftinvert_cpu.cpp
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
 */

#include "solve_eigval_cpu.h"
#include "solve_eigval_shiftinvert_op_general.h"

#include <cmath>

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_general_shiftinvert_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts)
{
    using OpB = Spectra::SparseSymMatProd<Precision>;

    constexpr Precision conditioning_tol = static_cast<Precision>(1e-10);
    constexpr int shift_steps = 50;

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = choose_ncv(n, k);

    logging::info(true, "Eigval (gen): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigvalPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    Precision shift = static_cast<Precision>(opts.sigma);
    ShiftInvertOpGeneral op(A, B, shift);
    OpB opB(B);

    // If the requested shift produces an almost singular shifted system, search
    // logarithmically for the smallest sufficiently stable negative shift.
    Precision conditioning = op.conditioning_ratio();
    if (conditioning < conditioning_tol) {
        logging::info(
            "Eigval (gen): ShiftInvert, sigma=", shift,
            " is too close to singular, adjusting..."
        );
        logging::down();

        for (int i = 0; i < shift_steps; ++i) {
            const Precision exponent = static_cast<Precision>(-12.0)
                                     + static_cast<Precision>(24.0 * i) / static_cast<Precision>(shift_steps - 1);

            shift = -std::pow(static_cast<Precision>(10), exponent);
            op.set_shift(shift);

            conditioning = op.conditioning_ratio();
            logging::info(
                true,
                "shift: ", shift,
                ", conditioning ratio: ", conditioning
            );

            if (conditioning >= conditioning_tol)
                break;
        }

        logging::up();
    }

    Spectra::SymGEigsShiftSolver<
        ShiftInvertOpGeneral,
        OpB,
        Spectra::GEigsMode::ShiftInvert
    > eig(op, opB, k, ncv_user, shift);

    eig.init();
    eig.compute(to_rule(opts.sort), opts.maxit, opts.tol);

    logging::error(eig.info() == Spectra::CompInfo::Successful,
                   "Spectra generalized (ShiftInvert) failed");

    DynamicVector evals = eig.eigenvalues();
    DynamicMatrix evecs = eig.eigenvectors();
    for (int i = 0; i < evals.size(); ++i)
        out.push_back({ static_cast<Precision>(evals[i]), evecs.col(i) });

    t.stop();
    logging::info(true, "Eigval (gen/shift-invert) elapsed: ", t.elapsed(), " ms");

    return out;
}
} // namespace fem::solver::detail
