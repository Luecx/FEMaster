/**
 * @file concentrated_heat_flux.h
 * @brief Defines concentrated nodal heat input for stationary conduction.
 *
 * Unlike surface heat density, CFLUX prescribes power directly at each target
 * node. The condition history and amplitude follow the common Condition model;
 * assembly adds power to the scalar thermal right-hand side without integration.
 */

#pragma once

#include "../condition.h"
#include "../../data/region.h"

namespace fem::bc {

/**
 * @brief Prescribes concentrated heat flow at every node of a compiled region.
 *
 * Each selected node receives the full nominal power, positive into the body.
 * The caller evaluates the shared amplitude and supplies its numerical scale.
 * This condition modifies only the scalar thermal RHS; the conductivity and
 * temperature constraints remain unchanged. Shared regions are retained by
 * ownership through the compiled NodeRegion pointer.
 */
struct ConcentratedHeatFlux : Condition {
    // Compiled nodes receiving the full power independently
    model::NodeRegion::Ptr region_ = nullptr;

    // Nominal heat flow per target node, in power units
    Precision power_ = Precision(0);

    // Add concentrated power to active scalar thermal equations
    void apply(
        model::ModelData&        model_data,
        model::Field&            load,
        TripletList*             load_tangent,
        constraint::Equations&   constraints,
        const model::SystemDofs& system_dof_ids,
        TripletList&             system_matrix,
        Precision                scale
    ) override;

    // Diagnostics
    std::string str() const override;
};

} // namespace fem::bc
