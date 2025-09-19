#pragma once
/**
 * @file cb_qr_rank.h
 * @brief Factor C_use with SparseQR and extract numerical rank & R block (in-place).
 *
 * Why in-place?
 * -------------
 * `Eigen::SparseQR` is non-copyable and non-movable. Returning a struct that
 * contains it by value and then assigning will fail to compile. Therefore this
 * API fills a caller-owned `QRRank` by reference.
 */

#include <Eigen/SparseQR>
#include "../../core/types_eig.h"

namespace fem { namespace constraint {

/// Holder for the QR factorization and rank info (non-copyable/non-movable).
struct QRRank {
    Eigen::SparseQR<SparseMatrix, Eigen::AMDOrdering<int>> sqr; ///< SparseQR factor of C_use
    SparseMatrix R;            ///< R = sqr.matrixR() (upper-triangular, may be rectangular)
    int          r{0};         ///< numerical rank
    Precision    R11_max_diag{0}; ///< max |diag(R11)| before thresholding

    QRRank() = default;
    QRRank(const QRRank&) = delete;
    QRRank& operator=(const QRRank&) = delete;
    QRRank(QRRank&&) = delete;
    QRRank& operator=(QRRank&&) = delete;
};

/**
 * @brief Factor C_use with SparseQR and compute numerical rank (in-place).
 *
 * @param C_use        m'×n_use matrix (after column compression).
 * @param rank_tol_rel Relative threshold for rank detection: pivot i is counted
 *                     while |R(i,i)| > rank_tol_rel * max|diag(R)|.
 * @param out          Output container filled with {sqr, R, r, R11_max_diag}.
 */
void factor_and_rank_inplace(const SparseMatrix& C_use,
                             Precision rank_tol_rel,
                             QRRank& out);

}} // namespace fem::constraint
