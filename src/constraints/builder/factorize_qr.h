/**
 * @file factorize_qr.h
 * @brief Declares wrappers around Eigen's sparse QR factorisation.
 *
 * The wrapper exposes configurable pivot thresholds and returns the factor in a
 * reusable form without moving the Eigen QR object.
 *
 * @see src/constraints/builder/factorize_qr.cpp
 * @see src/constraints/builder/preprocess.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>

namespace fem {
namespace constraint {

/**
 * @struct QRSettings
 * @brief Configures the sparse QR factorisation.
 */
struct QRSettings {
    Precision pivot_rel = 0.0; ///< Relative pivot threshold passed to SparseQR.
};

/**
 * @struct QRResult
 * @brief Holds the sparse QR factor and associated data.
 */
struct QRResult {
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr; ///< Factorisation handle.
    SparseMatrix R; ///< Upper-triangular factor `R`.
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P; ///< Column permutation.
};

/**
 * @brief Factorises the constraint matrix using sparse QR.
 *
 * @param C_use Constraint matrix to factorise.
 * @param settings QR configuration options.
 * @param out Result container populated in-place.
 */
void factorize_sparse_qr(const SparseMatrix& C_use, const QRSettings& settings, QRResult& out);

} // namespace constraint
} // namespace fem
