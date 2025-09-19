/******************************************************************************
* @file constraint_builder.h
 * @brief Build a reduced system map from linear constraints using QR.
 *
 * See .cpp for detailed algorithm notes.
 * @date    14.09.2025
 * @author  Finn Eggers
 ******************************************************************************/

#pragma once

#include "../constraint_set.h"
#include "../../core/types_eig.h"
#include <utility>
#include <string>
#include <vector>

namespace fem::constraint {

class ConstraintMap;

struct ConstraintBuilder {

    struct Options {
        Precision rank_tol_rel   = 1e-12;  ///< relative tol on diag(R) for rank
        Precision feas_tol_rel   = 1e-10;  ///< ||C u_p - d|| <= feas_tol_rel * ||d||
        bool      threshold_X    = true;   ///< currently unused (kept for API compat)
        Precision X_drop_tol_rel = 1e-12;  ///< currently unused (kept for API compat)
        int       suspect_rows_k = 10;     ///< currently unused (kept for API compat)
    };

    struct Report {
        Index m = 0;
        Index n = 0;
        Index rank = 0;
        bool  homogeneous = true;
        bool  feasible    = true;

        Precision R11_max_diag = 0;
        Precision residual_norm = 0;
        Precision d_norm        = 0;
        Index     n_redundant_rows = 0;

        std::vector<Index> slave_idx;    ///< size = rank
        std::vector<Index> master_idx;   ///< size = n - rank

        std::vector<Index> suspect_row_ids; ///< reserved
        std::string        log;             ///< reserved
    };

    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set);
    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set, const Options& opt);
};

} // namespace fem::constraint
