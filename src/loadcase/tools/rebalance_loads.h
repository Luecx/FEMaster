//
// Created by f_eggers on 10.02.2026.
//

#ifndef COMPILE_SH_RIGID_BODY_REBALANCING_H
#define COMPILE_SH_RIGID_BODY_REBALANCING_H

#include "../../model/model_data.h"

namespace fem {

/**
 * @brief Rebalances a nodal load field such that global resultant force and moment
 *        about the center of gravity vanish (free-free consistency), WITHOUT
 *        introducing constraints and WITHOUT recomputing accelerations (unlike
 *        inertia relief).
 *
 * This is a pure RHS correction limited to the rigid-body subspace (6 DOF).
 */
void rebalance_loads(model::ModelData& model_data, model::Field& global_load_mat);

} // namespace fem

#endif // COMPILE_SH_RIGID_BODY_REBALANCING_H