/******************************************************************************
* @file constraint_set.h
 * @brief Assemble the sparse constraint system `C u = d` from high-level equations.
 *
 * Scaling is disabled by default; we assemble exactly what is specified.
 * @date    14.09.2025                                                    @author Finn
 ******************************************************************************/

#pragma once

#include "../core/types_cls.h"
#include "../core/types_eig.h"
#include "equation.h"

#include <vector>

namespace fem::constraint {

struct ConstraintSet {
    // Input
    Equations           equations;
    const SystemDofIds* dof_map = nullptr;

    struct Options {
        bool      scale_columns     = false; ///< default off per your preference
        bool      scale_rows        = false; ///< default off per your preference
        Precision zero_row_drop_tol = 0;
    } opt;

    // Output
    Index          m = 0;
    Index          n = 0;
    SparseMatrix   C;
    DynamicVector  d;

    std::vector<Index> kept_row_ids;
    DynamicVector  col_scale;
    DynamicVector  row_scale;

    void assemble(const SystemDofIds& dofs, Index n_dofs);
};

} // namespace fem::constraint
