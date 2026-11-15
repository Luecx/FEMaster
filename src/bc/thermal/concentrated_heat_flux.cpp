/**
 * @file concentrated_heat_flux.cpp
 * @brief Assembles prescribed nodal power into the thermal balance.
 *
 * Concentrated heat flow acts directly on temperature DOFs and requires no
 * surface or volume Jacobian. The thermal analysis supplies amplitude scaling;
 * this implementation validates active targets and accumulates their power.
 */

#include "concentrated_heat_flux.h"

#include "../../core/logging.h"
#include "../../data/field.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Adds scaled power to every selected node in the scalar thermal RHS.
 *
 * Positive power supplies heat, negative power removes it. Nodes without a
 * thermal equation are rejected so prescribed power cannot be silently lost.
 * Constraint equations and sparse matrix contributions remain unchanged.
 *
 * @param load Nodal scalar heat-flow field modified in place.
 * @param system_dof_ids Scalar thermal numbering used to validate targets.
 * @param scale Amplitude or step multiplier evaluated by the analysis.
 */
void ConcentratedHeatFlux::apply(
    model::ModelData&,
    model::Field&          load,
    TripletList*,
    constraint::Equations&,
    const SystemDofIds&    system_dof_ids,
    TripletList&,
    Precision              scale
) {
    // Validate the thermal field and effective power before accumulating loads
    const Precision power = scale * power_;
    logging::error(region_ != nullptr,
        "CFLUX: target node region is not set");
    logging::error(load.domain == model::FieldDomain::NODE && load.components == 1,
        "CFLUX: target field must be a scalar NODE field");
    logging::error(static_cast<Index>(system_dof_ids.rows()) == load.rows && system_dof_ids.cols() == 1,
        "CFLUX: thermal DOF numbering must match the scalar nodal domain");
    logging::error(std::isfinite(power_) && std::isfinite(power),
        "CFLUX: prescribed and scaled power must be finite");

    // Concentrated power is already a nodal load and needs no geometric measure
    for (const ID node : region_->data()) {
        logging::error(node >= 0 && static_cast<Index>(node) < load.rows,
            "CFLUX: node is outside the thermal RHS domain: ", node);
        logging::error(system_dof_ids(node, 0) >= 0,
            "CFLUX: node has no active temperature DOF: ", node);
        load(node, 0) += power;
    }
}

std::string ConcentratedHeatFlux::str() const {
    std::ostringstream os;
    os << "CFLUX: target=" << (region_ ? region_->name : "?") << ", power=" << power_;
    return os.str();
}

} // namespace fem::bc
