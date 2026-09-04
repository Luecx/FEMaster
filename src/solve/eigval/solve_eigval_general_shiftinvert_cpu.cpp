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
#include <array>
#include <cmath>

namespace fem::solver::detail {

static Precision estimate_lambda_scale(
    const SparseMatrix& A,
    const SparseMatrix& B)
{
    std::array<Precision, 3> estimates{};
    int count = 0;

    for (int j = 0; j < 3; ++j) {
        const DynamicVector q  = DynamicVector::Random(A.rows());
        const DynamicVector Aq = A * q;
        const DynamicVector Bq = B * q;

        const Precision num = q.dot(Aq);
        const Precision den = q.dot(Bq);

        if (std::isfinite(num) &&
            std::isfinite(den) &&
            num > Precision(0) &&
            den > Precision(0))
        {
            estimates[count++] = num / den;
        }
    }

    if (count == 0)
        return Precision(0);

    std::sort(estimates.begin(), estimates.begin() + count);
    return estimates[count / 2];
}

std::vector<EigvalPair>
eigval_general_shiftinvert_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts)
{
    using OpB = Spectra::SparseSymMatProd<Precision>;

    constexpr Precision conditioning_tol = static_cast<Precision>(1e-10);

    const int n        = static_cast<int>(A.rows());
    const int ncv_user = choose_ncv(n, k);

    logging::info(true, "Eigval (gen): ShiftInvert, sigma=", opts.sigma, ", ncv=", ncv_user);

    std::vector<EigvalPair> out;
    out.reserve(static_cast<size_t>(k));

    Timer t{}; t.start();

    Precision shift = static_cast<Precision>(opts.sigma);
    ShiftInvertOpGeneral op(A, B, shift);
    OpB opB(B);

    // If the requested shift produces an almost singular shifted system, move
    // sigma to the negative side of the spectrum. The scale is estimated from
    // generalized Rayleigh quotients and increased geometrically if necessary.
    Precision conditioning = op.conditioning_ratio();
    if (conditioning < conditioning_tol) {
        logging::info(
            "Eigval (gen): ShiftInvert, sigma=", shift,
            " is too close to singular, adjusting..."
        );

        const Precision lambda_scale = std::abs(estimate_lambda_scale(A, B));
        Precision shift_factor = static_cast<Precision>(1e-4);

        logging::down();
        logging::info("estimated eigenvalue scale = ", lambda_scale);

        for (int i = 0; i < 8; ++i) {
            shift = -lambda_scale * shift_factor;
            op.set_shift(shift);

            conditioning = op.conditioning_ratio();
            logging::info(
                true,
                "shift: ", shift,
                ", conditioning ratio: ", conditioning
            );

            if (conditioning >= conditioning_tol)
                break;

            shift_factor *= static_cast<Precision>(10);
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
