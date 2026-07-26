/**
 * @file rebalance_loads.h
 * @brief Declares rigid-body rebalancing of assembled nodal loads.
 *
 * The rigid-body rebalancing utility modifies an assembled nodal load field
 * such that its resultant force and resultant moment vanish.
 *
 * The correction is restricted to six global rigid-load patterns:
 *
 * - three translational force patterns,
 * - three rotational force patterns.
 *
 * In contrast to inertia relief, the routine does not determine rigid-body
 * accelerations and does not construct inertia forces from a mass matrix.
 * Instead, it directly corrects the assembled right-hand side while leaving
 * the model constraints, stiffness matrix and active degrees of freedom
 * unchanged.
 *
 * The center of gravity required for the moment balance is obtained from the
 * structural elements of the model using the same integration convention as
 * the structural inertia-relief implementation.
 *
 * @see model::ModelData
 * @see model::Field
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#pragma once

#include "../../model/model_data.h"

namespace fem {

// Rigid-body load rebalancing. The assembled load field is modified in place
// by adding nodal forces whose total force and total moment cancel the original
// global resultants.
//
// Existing nodal moments participate in the resultant calculation when the
// field contains at least six components. The generated correction itself is
// applied only to the three translational load components.
void rebalance_loads(model::ModelData& model_data, model::Field& global_load_mat);

} // namespace fem