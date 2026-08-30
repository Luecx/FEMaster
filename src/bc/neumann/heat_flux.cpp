/**
 * @file heat_flux.cpp
 * @brief Implements prescribed surface heat-flux boundary conditions.
 *
 * Prescribed heat flux is a pure thermal right-hand-side contribution. The
 * scalar surface integration path multiplies the flux by the surface shape
 * functions and physical area measure before scattering it into the
 * one-component nodal thermal field.
 *
 * @see heat_flux.h
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "heat_flux.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the prescribed heat flux over the selected boundary surfaces.
 *
 * For constant heat flux q, each surface contributes the consistent nodal source
 * vector
 *
 *     f_q = integral_Gamma_q N^T q dGamma.
 *
 * Positive values are defined as heat entering the model. The optional LHS
 * assembly objects are unused because prescribed heat flux is independent of
 * temperature.
 *
 * @param model_data Compiled surface topology and reference nodal positions.
 * @param rhs One-component nodal thermal right-hand side.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 * @param system_dof_ids Unused optional system DOF map.
 * @param lhs Unused optional system-matrix triplet list.
 */
void HeatFlux::apply(model::ModelData&       model_data,
                     model::Field&           rhs,
                     Precision               time,
                     bool                    ignore_amplitude,
                     const SystemDofIds*      system_dof_ids,
                     TripletList*             lhs) {
    (void) system_dof_ids;
    (void) lhs;

    logging::error(region_ != nullptr,
        "HEATFLUX: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HEATFLUX: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "HEATFLUX: target field must be a NODE field with exactly one component");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision value = heat_flux_ * scale;

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "HEATFLUX: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "HEATFLUX: surface ", surface_id, " is not initialized");

        // Distribute the prescribed scalar flux with the surface interpolation
        // and complete physical area measure.
        surface->integrate_scalar_field(
            *model_data.positions_reference,
            rhs,
            [value](const Vec3&) -> Precision { return value; }
        );
    }
}

/**
 * Builds the diagnostic representation of the prescribed heat flux.
 *
 * @return Human-readable load description.
 */
std::string HeatFlux::str() const {
    std::ostringstream os;
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", value=" << heat_flux_;

    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
