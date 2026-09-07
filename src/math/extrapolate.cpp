/**
 * @file extrapolate.cpp
 * @brief Implements polynomial extrapolation between points in natural coordinates.
 *
 * @author Finn Eggers
 * @date 07.09.2026
 */

#include "extrapolate.h"

#include "../core/logging.h"

#include <Eigen/Cholesky>

namespace fem::math {

/**
 * Builds a linear operator that reconstructs a polynomial field from source
 * values and evaluates that field at target points.
 *
 * For the selected basis functions `F_j`, the source matrix is
 *
 *     A_ij = F_j(r_i, s_i, t_i).
 *
 * The polynomial coefficients are obtained by the normal equations
 *
 *     (A^T A) a = A^T q_source.
 *
 * Evaluating the same basis at the target points gives a matrix `B`, so the
 * complete source-to-target operator is
 *
 *     E = B (A^T A)^-1 A^T,
 *
 * and `q_target = E q_source`. All coordinates are natural/reference-space
 * coordinates; physical element geometry does not enter the recovery operator.
 *
 * @param source_points Points carrying the known values, one `(r,s,t)` per row.
 * @param target_points Points where values shall be reconstructed.
 * @param basis Polynomial basis defining the reconstruction assumption.
 * @return Source-to-target extrapolation matrix.
 */
RowMatrix extrapolate(const RowMatrix&                          source_points,
                      const RowMatrix&                          target_points,
                      std::initializer_list<ExtrapolationBasis> basis) {
    logging::error(source_points.cols() >= 3,
        "Extrapolation source points require at least three coordinates");
    logging::error(target_points.cols() >= 3,
        "Extrapolation target points require at least three coordinates");
    logging::error(basis.size() > 0,
        "Extrapolation requires at least one basis function");
    logging::error(source_points.rows() >= static_cast<Eigen::Index>(basis.size()),
        "Extrapolation requires at least as many source points as basis functions (",
        source_points.rows(), " vs ", basis.size(), ")");

    // Evaluate one selected monomial at one natural-coordinate point.
    const auto evaluate = [](ExtrapolationBasis function, Precision r, Precision s, Precision t) {
        switch (function) {
            case ExtrapolationBasis::F1:    return Precision(1);

            case ExtrapolationBasis::FR:    return r;
            case ExtrapolationBasis::FS:    return s;
            case ExtrapolationBasis::FT:    return t;

            case ExtrapolationBasis::FRR:   return r * r;
            case ExtrapolationBasis::FSS:   return s * s;
            case ExtrapolationBasis::FTT:   return t * t;

            case ExtrapolationBasis::FRS:   return r * s;
            case ExtrapolationBasis::FRT:   return r * t;
            case ExtrapolationBasis::FST:   return s * t;
            case ExtrapolationBasis::FRST:  return r * s * t;

            case ExtrapolationBasis::FRRS:  return r * r * s;
            case ExtrapolationBasis::FRRT:  return r * r * t;
            case ExtrapolationBasis::FSSR:  return s * s * r;
            case ExtrapolationBasis::FSST:  return s * s * t;
            case ExtrapolationBasis::FTTR:  return t * t * r;
            case ExtrapolationBasis::FTTS:  return t * t * s;

            case ExtrapolationBasis::FRRST: return r * r * s * t;
            case ExtrapolationBasis::FRSST: return r * s * s * t;
            case ExtrapolationBasis::FRSTT: return r * s * t * t;
        }

        return Precision(0);
    };

    const Eigen::Index n_source = source_points.rows();
    const Eigen::Index n_target = target_points.rows();
    const Eigen::Index n_basis  = static_cast<Eigen::Index>(basis.size());

    RowMatrix A(n_source, n_basis);
    RowMatrix B(n_target, n_basis);

    // Assemble the polynomial evaluation matrices at source and target points.
    Eigen::Index j = 0;
    for (const auto function : basis) {
        for (Eigen::Index i = 0; i < n_source; ++i)
            A(i, j) = evaluate(function, source_points(i, 0), source_points(i, 1), source_points(i, 2));

        for (Eigen::Index i = 0; i < n_target; ++i)
            B(i, j) = evaluate(function, target_points(i, 0), target_points(i, 1), target_points(i, 2));

        ++j;
    }

    // Solve the normal system once for the complete source-to-coefficient map.
    const DynamicMatrix At     = A.transpose();
    const DynamicMatrix normal = At * A;

    Eigen::LDLT<DynamicMatrix> solver(normal);
    logging::error(solver.info() == Eigen::Success,
        "Extrapolation normal-system factorization failed");

    const DynamicMatrix reconstruction = solver.solve(At);
    logging::error(solver.info() == Eigen::Success && reconstruction.allFinite(),
        "Extrapolation normal-system solution failed");

    return RowMatrix(B * reconstruction);
}

} // namespace fem::math
