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

#include <algorithm>
#include <cmath>

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_general_shiftinvert_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts)
{
    using OpB = Spectra::SparseSymMatProd<Precision>;

    constexpr Precision conditioning_tol = static_cast<Precision>(1e-11);

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = choose_ncv(n, k);

    logging::info(true, "Eigval (gen): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigvalPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    Precision shift = static_cast<Precision>(opts.sigma);
    ShiftInvertOpGeneral op(A, B, shift);
    OpB opB(B);

    // If the requested shift is nearly singular, use the externally supplied
    // characteristic eigenvalue scale to choose a small negative starting shift.
    Precision conditioning = op.conditioning_ratio();
    if (conditioning < conditioning_tol) {
        logging::info(
            "Eigval (gen): ShiftInvert, sigma=", shift,
            " is too close to singular, adjusting..."
        );
        logging::up();

        const Precision lambda_scale = std::abs(opts.lambda_scale);
        const bool has_scale = std::isfinite(lambda_scale) && lambda_scale > Precision(0);

        Precision shift_magnitude = has_scale
            ? std::max(lambda_scale * static_cast<Precision>(1e-4), static_cast<Precision>(1e-12))
            : static_cast<Precision>(1e-12);

        const Precision growth = static_cast<Precision>(100);
        const int max_attempts = has_scale ? 8 : 13;

        if (has_scale)
            logging::info(true, "Rayleigh eigenvalue scale: ", lambda_scale);

        for (int i = 0; i < max_attempts; ++i) {
            shift = -shift_magnitude;
            op.set_shift(shift);

            conditioning = op.conditioning_ratio();
            logging::info(
                true,
                "shift: ", shift,
                ", conditioning ratio: ", conditioning
            );

            if (conditioning >= conditioning_tol)
                break;

            shift_magnitude *= growth;
        }

        logging::down();
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
