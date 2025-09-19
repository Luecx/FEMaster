/******************************************************************************
* @file build_x.h
 * @brief Build sparse X = -R11^{-1}R12 as column lists via upper-tri backsolves.
 ******************************************************************************/
#pragma once
#include "../../core/types_eig.h"
#include <vector>

namespace fem { namespace constraint {

struct XCols {
    // For each master j (local 0..nm_use-1), store list of (i, val) for X(i,j)
    std::vector<std::vector<std::pair<int, Precision>>> cols;
};

XCols build_X_cols_from_R(const SparseMatrix& R, int r);

}} // namespace
