/**
 * @file reduce_mat_to_vec.h
 * @brief Converts between nodal fields and active system vectors.
 *
 * The helpers operate on arbitrary node-by-component `SystemDofIds` mappings.
 * Negative entries denote inactive components; non-negative entries address the
 * corresponding coefficient of the compact algebraic system vector.
 *
 * @author Created by Finn Eggers (c)
 * all rights reserved
 * @date Created on 28.08.2024
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/field.h"

namespace fem { namespace mattools {

/**
 * @brief Reduces a nodal field into the active algebraic system vector.
 *
 * @param dof_ids Node-by-component global system identifiers.
 * @param field Nodal field containing the corresponding physical components.
 * @return Compact vector addressed by the non-negative entries of `dof_ids`.
 */
DynamicVector reduce_mat_to_vec(
    const SystemDofIds& dof_ids,
    const model::Field& field
);

/**
 * @brief Expands an active algebraic system vector into nodal field storage.
 *
 * @param dof_ids Node-by-component global system identifiers.
 * @param reduced_vector Compact vector containing active system coefficients.
 * @return NODE-domain field with the same component count as `dof_ids`.
 */
model::Field expand_vec_to_mat(
    const SystemDofIds& dof_ids,
    const DynamicVector& reduced_vector
);

} } // namespace fem::mattools
