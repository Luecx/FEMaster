/******************************************************************************
* @file factorize_qr.h
* @brief SparseQR factorization wrapper with configurable pivot threshold.
******************************************************************************/
#pragma once
#include "../../core/types_eig.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

namespace fem { namespace constraint {

struct QRSettings {
    Precision pivot_rel = 0.0;  ///< 0..0.1 typical; interpreted by SparseQR
};

struct QRResult {
    // Keep the handle so callers can call qr.solve later.
    Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr; ///< handle
    SparseMatrix R;    ///< matrixR() copy (upper-triangular)
    Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> P; ///< cols permutation
};

/**
 * @brief Factorize C_use into QR; writes into @p out to avoid moving SparseQR.
 */
void factorize_sparse_qr(const SparseMatrix& C_use, const QRSettings& s, QRResult& out);

}} // namespace fem::constraint
