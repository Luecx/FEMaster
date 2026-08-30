/**
 * @file load_p.cpp
 * @brief Implements pressure integration along geometric surface normals.
 *
 * Each selected surface converts the scalar pressure into a vector traction at
 * its quadrature points using the physical surface normal and distributes the
 * integrated contribution consistently to its nodes. Optional amplitude scaling
 * is evaluated once before the surface loop.
 *
 * @see load_p.h
 * @see SurfaceInterface
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_p.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

/**
 * Integrates scalar pressure over the selected surfaces.
 *
 * The physical surface normal converts the scalar pressure into a global
 * traction vector. Amplitude scaling is evaluated once, while surface
 * integration performs physical area weighting and consistent nodal
 * distribution. Positive pressure acts opposite to the surface normal.
 *
 * @param model_data Global nodal positions required for surface geometry.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void PLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    // Validate the geometry, selected region and structural target field before
    // evaluating pressure magnitude or traversing any surfaces.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "PLoad: target surface region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "PLoad: target field must be a NODE field with at least six components");

    // Pressure is scalar, so its temporal multiplier can be evaluated once and
    // folded into the nominal magnitude before geometric integration.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision pressure = pressure_ * scale;
    const auto& positions    = *model_data.positions;

    // Integrate pressure independently over every selected compiled surface.
    for (const ID surface_id : *region_) {
        // Empty compiled surface slots do not contribute to the load field.
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        // The surface performs quadrature and nodal shape-function weighting.
        // The callback supplies the pressure traction at the current global
        // integration-point position.
        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                // Recover the natural coordinates needed by the normal
                // implementation and apply positive pressure opposite to the
                // oriented physical surface normal.
                const Vec2 local = surface->global_to_local(position, positions);
                return -pressure * surface->normal(positions, local);
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the pressure load.
 *
 * The result identifies the target surface region, its size and the nominal
 * scalar pressure before temporal amplitude scaling.
 *
 * @return Human-readable load description.
 */
std::string PLoad::str() const {
    std::ostringstream os;

    // Report the surface set, its current size and the unscaled pressure.
    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", p=" << pressure_;

    // Include the amplitude only when temporal scaling is configured.
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
