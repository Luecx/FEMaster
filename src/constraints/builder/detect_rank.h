/******************************************************************************
* @file detect_rank.h
 * @brief Rank detection from diagonal of R with relative tolerance.
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"

namespace fem { namespace constraint {

struct RankSettings {
    Precision rel_tol = 1e-12;
};

struct RankInfo {
    int       r = 0;
    Precision max_abs_diag = 0;
};

RankInfo detect_rank_from_R(const SparseMatrix& R, const RankSettings& s);

}} // namespace
