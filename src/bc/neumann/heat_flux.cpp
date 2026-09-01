/**
 * @file heat_flux.cpp
 * @brief Implements consistent prescribed thermal surface heat-flux assembly.
 *
 * The implementation validates the scalar thermal RHS and selected compiled
 * surfaces, evaluates optional amplitude scaling once, and delegates reference-
 * surface quadrature and shape-function distribution to `SurfaceInterface`.
 * Each selected surface therefore contributes a consistent nodal heat-flow
 * vector rather than an element-nodal result field.
 *
 * @see heat_flux.h
 * @see model::SurfaceInterface
 *
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
 * Integrates the prescribed scalar heat flux over all selected surfaces.
 *
 * For each surface, the scalar callback is integrated with the surface shape
 * functions and physical reference-area measure. The resulting nodal heat-flow
 * contribution is accumulated directly in the supplied one-component NODE
 * field. The nominal flux is amplitude-scaled before surface traversal because
 * the same scalar multiplier applies to the complete region.
 *
 * @param model_data Compiled surfaces and reference nodal geometry.
 * @param rhs Scalar nodal thermal right-hand-side field receiving the heat flow.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether the nominal unscaled flux should be applied.
 */
void HeatFlux::apply(model::ModelData& model_data,
                     model::Field&     rhs,
                     Precision         time,
                     bool              ignore_amplitude) {
    // Validate the semantic target, reference geometry, thermal target layout and
    // prescribed physical value before entering surface integration.
    logging::error(region_ != nullptr,
        "HEATFLUX: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HEATFLUX: reference positions are not initialized");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components == 1,
        "HEATFLUX: target field must be a NODE field with exactly one component");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    // Evaluate the temporal amplitude once and fold it into the nominal scalar
    // flux before traversing the target surfaces.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision value = heat_flux_ * scale;

    // Integrate the same scalar flux over every compiled surface in the selected
    // region using the reference configuration as the surface measure.
    for (ID surface_id : *region_) {
        // Validate the dense compiled surface identifier before resolving the
        // corresponding geometry object.
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "HEATFLUX: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "HEATFLUX: surface ", surface_id, " is not initialized");

        // The surface integrator evaluates shape functions, reference-area
        // Jacobians and quadrature weights and scatters the consistent scalar
        // nodal contribution directly into `rhs`.
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
 * @return Human-readable target surface region, nominal flux and amplitude.
 */
std::string HeatFlux::str() const {
    std::ostringstream os;

    // Report the selected surface set, its current size and the unscaled nominal
    // heat-flux value.
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", value=" << heat_flux_;

    // Include temporal scaling only when an amplitude is assigned.
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
