/**
 * @file eigen_common.h
 * @brief Shared utilities, includes, and helpers for eigen solvers.
 *
 * @details Provides:
 *          - Spectra include set used by all solver variants.
 *          - Small helpers to choose subspace size (ncv), map sorting rules,
 *            and ensure Eigen sparse matrices are in compressed (CSC) format.
 *
 *          These functions are internal to the solver implementation and are
 *          not part of the public API.
 *
 * @author
 *   Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *   All rights reserved.
 * @date   Created on 19.09.2025
 */
#pragma once
#include "eigen.h"                // public API
#include "../core/logging.h"
#include "../core/timer.h"

#include <Eigen/SparseCore>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>

#include <Spectra/MatOp/SparseSymMatProd.h>

#include <algorithm>
#include <vector>
#include <cmath>

namespace fem::solver::detail {

/**
 * @brief Map high-level sorting rule to Spectra’s SortRule.
 */
inline Spectra::SortRule to_rule(EigenOpts::Sort s) {
    using SR   = Spectra::SortRule;
    using Sort = EigenOpts::Sort;
    switch (s) {
        case Sort::LargestAlge:  return SR::LargestAlge;
        case Sort::LargestMagn:  return SR::LargestMagn;
        case Sort::SmallestAlge: return SR::SmallestAlge;
        case Sort::SmallestMagn: default: return SR::SmallestMagn;
    }
}

/**
 * @brief Choose subspace size ncv heuristically: k < ncv ≤ n.
 */
inline int choose_ncv(int n, int k) {
    const int lo  = k + 2;
    const int hi  = 2 * k + 20;
    int ncv       = std::min(n, std::max(lo, hi));
    if (ncv <= k) ncv = std::min(n, k + 2);
    return ncv;
}

/**
 * @brief Ensure a sparse matrix is in compressed (CSC) storage.
 */
inline void ensure_compressed(const SparseMatrix& M_const) {
    auto& M = const_cast<SparseMatrix&>(M_const);
    if (!M.isCompressed()) M.makeCompressed();
}

} // namespace fem::solver::detail
