/**
 * @file builder.h
 * @brief Declares the constraint builder that produces null-space maps.
 *
 * The builder consumes an assembled constraint set and computes the null-space
 * representation along with diagnostic information.
 *
 * @see src/constraints/builder/builder.cpp
 * @see src/constraints/constraint_map.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "../../core/types_eig.h"
#include "../constraint_set.h"

#include <string>
#include <utility>
#include <vector>

namespace fem {
namespace constraint {

class ConstraintMap;

/**
 * @struct ConstraintBuilder
 * @brief Provides static helpers to construct constraint maps.
 */
struct ConstraintBuilder {
    /**
     * @struct Options
     * @brief Configuration parameters for the builder.
     */
    struct Options {
        Precision rank_tol_rel = 1e-12; ///< Relative tolerance on diag(R) for rank detection.
        Precision feas_tol_rel = 1e-10; ///< Relative tolerance for feasibility check `||C u_p - d||`.
        bool threshold_X = true;        ///< Reserved flag for optional thresholding of `X`.
        Precision X_drop_tol_rel = 1e-12; ///< Reserved tolerance for dropping small entries in `X`.
        int suspect_rows_k = 10;          ///< Reserved parameter for suspect-row reporting.
    };

    /**
     * @struct Report
     * @brief Diagnostics describing the builder outcome.
     */
    struct Report {
        Index m = 0;          ///< Number of constraint equations.
        Index n = 0;          ///< Number of DOFs in the original system.
        Index rank = 0;       ///< Numerical rank of the constraint system.
        bool homogeneous = true; ///< Whether `d` was (numerically) zero.
        bool feasible = true;    ///< Whether the system passed feasibility checks.

        Precision R11_max_diag = 0; ///< Maximum diagonal entry in the R11 block.
        Precision residual_norm = 0; ///< Norm of the residual `||C u_p - d||`.
        Precision d_norm = 0;        ///< Norm of the right-hand side `||d||`.
        Index n_redundant_rows = 0;  ///< Number of detected redundant rows.

        std::vector<Index> slave_idx;  ///< Indices mapped to slave DOFs (size `rank`).
        std::vector<Index> master_idx; ///< Indices mapped to master DOFs (size `n - rank`).

        std::vector<Index> suspect_row_ids; ///< Reserved list of suspect row indices.
        std::string log;                   ///< Reserved log output.
    };

    /**
     * @brief Builds a constraint map using default options.
     *
     * @param set Assembled constraint set to process.
     * @return std::pair<ConstraintMap, Report> Map and diagnostics.
     */
    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set);

    /**
     * @brief Builds a constraint map using the provided options.
     *
     * @param set Assembled constraint set to process.
     * @param opt Builder configuration parameters.
     * @return std::pair<ConstraintMap, Report> Map and diagnostics.
     */
    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set, const Options& opt);
};

} // namespace constraint
} // namespace fem
