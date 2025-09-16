/******************************************************************************
 * @file constraint_builder.h
 * @brief Build a reduced system map from linear constraints using QR.
 *
 * -----------------------------------------------------------------------------
 * ## What problem does this solve?
 *
 * In finite element analysis we often need to enforce linear constraints:
 *
 *      C u = d
 *
 * where:
 *   - `u ∈ R^n` are all displacement DOFs,
 *   - `C ∈ R^{m×n}` encodes the constraints (one row = one equation),
 *   - `d ∈ R^m` is the right-hand side (usually 0 for homogeneous).
 *
 * Examples:
 *   - Prescribing u0 = 0 (a support),
 *   - Forcing u3 - 2 u7 = 0 (a multipoint constraint),
 *   - Linking several nodes so their displacements add to a fixed value.
 *
 * -----------------------------------------------------------------------------
 * ## Why do we reduce?
 *
 * Solvers typically want an **unconstrained system** A q = b.
 * Directly solving with constraints is awkward because `C u = d` must hold
 * *exactly*. One way is to build a *map*:
 *
 *      u = u_p + T q
 *
 * where:
 *   - `q ∈ R^(n−r)` are the **independent (master) DOFs**,
 *   - the other DOFs (slaves) are expressed in terms of q,
 *   - `u_p` is a *particular solution* (needed if d ≠ 0).
 *
 * Plugging this into Ku=f gives:
 *
 *      Tᵀ K T q = Tᵀ (f − K u_p)
 *
 * which is a smaller, unconstrained system that can be solved normally.
 *
 * -----------------------------------------------------------------------------
 * ## How does ConstraintBuilder find T and u_p?
 *
 * 1. Start from C u = d.
 * 2. Run a rank-revealing QR factorization with column pivoting:
 *
 *        C P = Q R
 *
 *    where P permutes the columns (DOFs), Q is orthogonal, and R is upper-triangular.
 *
 * 3. From R we can separate:
 *      - the first `r = rank(C)` columns (slave DOFs),
 *      - the remaining `n−r` columns (master DOFs).
 *
 * 4. Solve for the slaves in terms of masters:
 *
 *        R11 u_S + R12 u_M = d'
 *        => u_S = −R11⁻¹ R12 u_M   +   (particular solution if d' ≠ 0)
 *
 *    This yields:
 *      - X = −R11⁻¹ R12, the linear relation between slave and master DOFs,
 *      - T = [I; X], the overall mapping from master DOFs q to full vector u,
 *      - u_p, the offset vector needed if d ≠ 0.
 *
 * -----------------------------------------------------------------------------
 * ## What do you get?
 *
 * After build():
 *   - `ConstraintMap` with:
 *       * `T`  (n×(n−r)) mapping masters → full vector,
 *       * `X`  (r×(n−r)) slave coefficients,
 *       * `u_p` (n) particular solution.
 *   - `Report` with diagnostics:
 *       * rank(C), number of masters/slaves,
 *       * whether constraints were feasible (consistent with d),
 *       * indices of suspect contradictory rows.
 *
 * -----------------------------------------------------------------------------
 * ## ELI5 summary:
 *
 * Think of constraints like puppet strings:
 *   - Some DOFs are **masters**: they move freely (like your hands).
 *   - Other DOFs are **slaves**: their motion is tied to the masters
 *     (like puppet legs tied to the strings).
 *   - If d ≠ 0, the puppet is shifted a bit before moving (that’s u_p).
 *
 * The builder figures out who is master/slave, how the strings are tied,
 * and whether the puppet can actually move the way the constraints demand.
 *
 * -----------------------------------------------------------------------------
 * @see constraint_builder.cpp
 * @see constraint_map.h / constraint_map.cpp
 * @date    14.09.2025
 * @author  Finn Eggers
 ******************************************************************************/


#pragma once

#include "constraint_set.h"
#include "../core/types_eig.h"
#include <utility>
#include <string>
#include <vector>

namespace fem::constraint {

/**
 * @brief Forward declaration of the target immutable map type.
 *
 * ConstraintMap stores the mapping u = u_p + T q and provides application/
 * assembly helpers. It is constructed exclusively via ConstraintBuilder.
 */
class ConstraintMap;

/**
 * @brief Builder that factors C and constructs a ConstraintMap and a Report.
 *
 * This class has no state; it exposes static build() entry points that return
 * a freshly constructed ConstraintMap together with a diagnostic Report.
 */
struct ConstraintBuilder {

    /**
     * @brief Numerical/assembly options controlling rank detection and sparsity.
     */
    struct Options {
        /// Relative rank tolerance against max |diag(R)|.
        /// A diagonal |R_kk| ≤ rank_tol_rel * max_diag is treated as zero.
        Precision rank_tol_rel   = 1e-12;

        /// Feasibility tolerance for inhomogeneous sets (used on ||C u_p − d||).
        /// The test is:  ||res|| ≤ feas_tol_rel * ||d||  (or *1 if ||d||=0).
        Precision feas_tol_rel   = 1e-10;

        /// If true, drop tiny entries in X columns relative to that column’s norm.
        bool      threshold_X    = true;

        /// Relative drop tolerance for X entries (per-column).
        Precision X_drop_tol_rel = 1e-12;

        /// Number of largest-residual rows to report when infeasible.
        int       suspect_rows_k = 10;
    };

    /**
     * @brief Diagnostics and artifacts from the build() process.
     *
     * This struct summarizes rank, feasibility, partitions, and a human-readable
     * log line. Index vectors reflect the final slave/master ordering (global
     * DOF indices).
     */
    struct Report {
        Index m = 0;                  ///< Number of constraint rows provided.
        Index n = 0;                  ///< Number of DOFs (columns of C).
        Index rank = 0;               ///< Numerical rank detected for C.
        bool  homogeneous = true;     ///< True if d is (numerically) zero.
        bool  feasible    = true;     ///< Feasible inhomogeneous system (lsq residual small).

        Precision R11_max_diag = 0;   ///< Max |diag(R)| used for relative rank tolerance.
        Precision residual_norm = 0;  ///< ||C u_p − d|| (0 for homogeneous).
        Precision d_norm        = 0;  ///< ||d|| (0 for homogeneous).
        Index     n_redundant_rows = 0; ///< Estimated m - rank when feasible.

        std::vector<Index> slave_idx;   ///< Size = rank. Global DOF ids treated as slaves.
        std::vector<Index> master_idx;  ///< Size = n - rank. Global DOF ids treated as masters.

        std::vector<Index> suspect_row_ids; ///< Top-K residual rows if infeasible (indices into original equations).
        std::string        log;             ///< Summary string suitable for logging.
    };

    /**
     * @brief Build a ConstraintMap with default options.
     *
     * @param set Assembled ConstraintSet (contains C, d, sizes, and bookkeeping).
     * @return pair {ConstraintMap, Report}
     */
    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set);

    /**
     * @brief Build a ConstraintMap with user-specified options.
     *
     * This overload allows tuning rank/feasibility tolerances and sparsity control.
     *
     * @param set Assembled ConstraintSet (contains C, d, sizes, and bookkeeping).
     * @param opt Numerical/assembly options (rank tolerance, X thresholding, etc.).
     * @return pair {ConstraintMap, Report}
     */
    static std::pair<class ConstraintMap, Report> build(const ConstraintSet& set, const Options& opt);
};

} // namespace fem::constraint
