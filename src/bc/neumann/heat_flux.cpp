/**
 * @file heat_flux.cpp
 * @brief Implements prescribed surface heat-flux boundary conditions.
 *
 * Prescribed heat flux is a pure thermal right-hand-side contribution. FEMaster
 * intentionally reuses the existing three-component surface-vector integration
 * path: the scalar heat-flow contribution is placed in component zero while the
 * remaining components stay zero. No thermal-specific surface integrator is
 * required.
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
 * For constant heat flux \f$q\f$, each surface contributes the consistent nodal
 * source vector \f$\int_{\Gamma_q}\mathbf{N}^T q\,\mathrm{d}\Gamma\f$. Positive
 * values are defined as heat entering the model. The optional LHS assembly
 * objects are unused because prescribed heat flux is independent of temperature.
 *
 * @param model_data Compiled surface topology and reference nodal positions.
 * @param rhs Three-component nodal thermal RHS. Only component zero is modified.
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
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 3,
        "HEATFLUX: target field must be a NODE field with at least three components");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Vec3 value{heat_flux_ * scale, Precision(0), Precision(0)};

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "HEATFLUX: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "HEATFLUX: surface ", surface_id, " is not initialized");

        surface->integrate_vector_field(
            *model_data.positions_reference,
            rhs,
            [value](const Vec3&) -> Vec3 { return value; }
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
