// ----------------------------------------------------------------------------
#pragma once

#include "../../core/types_eig.h"
#include <vector>
#include <utility>

namespace fem { namespace constraint {

struct SimpleRow {
    int       row;   // row index in original C
    int       col;   // the only column with nonzero
    Precision a;     // C(row,col)
    Precision d;     // set.d[row] (or 0 if d is empty)
};

struct RowEntry { int k; Precision a; }; // for R11 row-wise storage

}} // namespace fem::constraint

// ============================================================================
