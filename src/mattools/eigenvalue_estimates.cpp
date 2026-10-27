/**
 * @file eigenvalue_estimates.cpp
 * @brief Implements directional and diagonal estimates of generalized eigenvalues.
 *
 * The Rayleigh quotient measures one vector in A x = lambda B x. Vector
 * normalization avoids dependence on its amplitude; a relative denominator
 * check handles cancellation when B is indefinite. The diagonal estimate
 * supplies a characteristic eigenvalue scale without solving an eigenproblem.
 *
 * Solver configuration, physical interpretation and fallback decisions belong
 * to the calling analysis. No factorization or model state is modified.
 *
 * @see eigenvalue_estimates.h
 * @author Finn Eggers
 * @date 27.10.2026
 */

#include "eigenvalue_estimates.h"

#include <cmath>
#include <limits>

namespace fem {
namespace mattools {

/**
 * Evaluates the generalized Rayleigh quotient in the supplied coordinates.
 *
 * Normalize x before matrix products so a rescaled trial vector describes the
 * same direction without squaring its amplitude. The denominator may be
 * negative for indefinite B. Reject it only when it is nonfinite or negligible
 * relative to ||B x||_2 for the normalized direction, which detects cancellation
 * without comparing quantities having different physical units.
 *
 * @return A finite quotient, including zero or negative values, or NaN when
 *         dimensions, vector data or arithmetic do not define a reliable quotient.
 */
Precision estimate_eigenvalue_rayleigh(const SparseMatrix& A, const SparseMatrix& B, const DynamicVector& x) {
    // Require a common nonempty coordinate space before any matrix product.
    const Precision invalid = std::numeric_limits<Precision>::quiet_NaN();
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()
        || x.size() != A.rows() || x.size() == 0 || !x.allFinite()) {
        return invalid;
    }

    // Use a unit direction; stableNorm avoids squaring very large/small entries.
    const Precision vector_norm = x.stableNorm();
    if (!std::isfinite(vector_norm) || vector_norm == Precision(0)) return invalid;

    const DynamicVector direction = x / vector_norm;
    const DynamicVector Ax        = A * direction;
    const DynamicVector Bx        = B * direction;
    const Precision numerator         = direction.dot(Ax);
    const Precision denominator       = direction.dot(Bx);
    const Precision denominator_scale = Bx.stableNorm();

    // An indefinite B can have x^T B x = 0 despite B x != 0. Treat such a
    // direction as undefined rather than reporting an arbitrarily large shift.
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon();
    if (!std::isfinite(numerator) || !std::isfinite(denominator) || !std::isfinite(denominator_scale)
        || std::abs(denominator) <= tolerance * denominator_scale) {
        return invalid;
    }

    const Precision quotient = numerator / denominator;
    return std::isfinite(quotient) ? quotient : invalid;
}

/**
 * Estimates a generalized eigenvalue scale from absolute diagonal sums.
 *
 * The ratio has the units of lambda in A x = lambda B x and scales linearly
 * with A and inversely with B. It does not measure the minimum eigenvalue:
 * positive diagonal sums can coexist with rigid-body or other null modes.
 * A zero numerator is preserved; the calling analysis decides whether it needs
 * a nonzero fallback. Only the diagonals are inspected, without changing storage.
 *
 * @return A finite nonnegative scale, or NaN for incompatible matrices,
 *         nonfinite diagonal sums or a zero denominator.
 */
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B) {
    // Validate dimensions before accessing the common eigenvalue coordinates.
    const Precision invalid = std::numeric_limits<Precision>::quiet_NaN();
    if (A.rows() != A.cols() || B.rows() != B.cols() || A.rows() != B.rows()) return invalid;

    // Absolute sums measure magnitude even when diagonal signs differ.
    const Precision numerator   = A.diagonal().cwiseAbs().sum();
    const Precision denominator = B.diagonal().cwiseAbs().sum();
    if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator == Precision(0)) return invalid;

    const Precision scale = numerator / denominator;
    return std::isfinite(scale) ? scale : invalid;
}

} // namespace mattools
} // namespace fem
