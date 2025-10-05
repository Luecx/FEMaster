/******************************************************************************
 * @file detect_rank.h
 * @brief Declares helpers that determine the numerical rank from QR factors.
 *
 * Rank detection inspects the diagonal of the upper-triangular factor using a
 * relative tolerance to distinguish significant pivots.
 *
 * @see src/constraints/builder/detect_rank.cpp
 * @see src/constraints/builder/factorize_qr.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "../../core/types_eig.h"

namespace fem {
namespace constraint {

/******************************************************************************
 * @struct RankSettings
 * @brief Parameters controlling rank detection.
 ******************************************************************************/
struct RankSettings {
    Precision rel_tol = 1e-12; ///< Relative tolerance applied to diagonal entries.
};

/******************************************************************************
 * @struct RankInfo
 * @brief Holds the detected rank and auxiliary metrics.
 ******************************************************************************/
struct RankInfo {
    int r = 0;                  ///< Detected rank.
    Precision max_abs_diag = 0; ///< Maximum absolute diagonal entry.
};

/******************************************************************************
 * @brief Detects the numerical rank from the upper-triangular factor `R`.
 *
 * @param R Upper-triangular factor (typically from sparse QR).
 * @param settings Rank detection settings.
 * @return RankInfo Detected rank and diagnostics.
 ******************************************************************************/
RankInfo detect_rank_from_R(const SparseMatrix& R, const RankSettings& settings);

} // namespace constraint
} // namespace fem
