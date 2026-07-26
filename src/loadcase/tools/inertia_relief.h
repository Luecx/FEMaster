/**
 * @file inertia_relief.h
 * @brief Declares inertia relief for unconstrained structural systems.
 *
 * Inertia relief balances the externally applied force and moment resultants
 * of a free or insufficiently constrained structural model by superimposing
 * equivalent rigid-body inertia loads.
 *
 * The required translational and angular accelerations are determined from the
 * total mass, center of gravity and inertia tensor of the model. Distributed
 * structural mass is always considered. Concentrated point masses and their
 * rotary inertias may optionally be included.
 *
 * The resulting inertia loads are assembled through the existing
 * `bc::InertialLoad` implementation and added directly to the supplied nodal
 * load field.
 *
 * The routine modifies only the assembled load field. It does not alter:
 *
 * - the stiffness matrix,
 * - the mass matrix,
 * - model constraints,
 * - active degrees of freedom,
 * - nodal positions,
 * - element state.
 *
 * @see bc::InertialLoad
 * @see feature::PointMass
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 26.07.2026
 */

#pragma once

#include "../../model/model_data.h"

namespace fem {

// Apply inertia relief to an assembled nodal load field.
//
// The routine determines the translational and angular rigid-body
// accelerations required to balance the current global force and moment
// resultants. The corresponding signed inertia loads are assembled through
// bc::InertialLoad and added directly to global_load_mat.
//
// Distributed structural mass is always included. Concentrated point masses,
// their offsets from the center of gravity and their rotary inertias are
// included when consider_point_masses is true.
void apply_inertia_relief(
    model::ModelData& model_data,
    model::Field&     global_load_mat,
    bool              consider_point_masses = true
);

} // namespace fem

