/**
 * @file reduce_mat_to_vec.h
 * @brief Declares conversion between nodal fields and active system vectors.
 *
 * The helpers operate on arbitrary node-by-component `SystemDofIds` mappings.
 * Negative entries denote inactive components; non-negative entries address the
 * corresponding coefficient of the compact algebraic system vector. The same
 * interface therefore supports six-component structural mappings and scalar
 * node-by-one thermal mappings without assuming a fixed field width.
 *
 * Detailed validation and mapping semantics are documented at the definitions in
 * `reduce_mat_to_vec.cpp`.
 *
 * @author Created by Finn Eggers (c)
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/field.h"

namespace fem { namespace mattools {

// Convert a NODE field into compact active-system ordering. Every non-negative
// entry in `dof_ids` identifies the exact destination coefficient; inactive
// entries are skipped rather than assumed to occur in any particular pattern.
DynamicVector reduce_mat_to_vec(
    const SystemDofIds& dof_ids,
    const model::Field& field
);

// Expand a compact active-system vector back into NODE storage with the exact row
// and component dimensions of `dof_ids`. Inactive components remain zero.
model::Field expand_vec_to_mat(
    const SystemDofIds& dof_ids,
    const DynamicVector& reduced_vector
);

} } // namespace fem::mattools
