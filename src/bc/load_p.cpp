/**
 * @file load_p.cpp
 * @brief Implements pressure integration along geometric surface normals.
 *
 * Each selected surface converts the scalar pressure into a vector traction at
 * its quadrature points, using its current geometric normal, and distributes
 * the integrated contribution consistently to its nodes. Optional amplitude
 * scaling is evaluated once before the surface loop.
 *
 * @see load_p.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_p.h"

#include "../core/logging.h"
#include "../model/model_data.h"

#include <sstream>

namespace fem {
namespace bc {

/**
 * Integrates scalar pressure over the selected surfaces.
 *
 * The current surface normal converts pressure into a global traction vector.
 * Amplitude scaling is evaluated once, while surface integration performs
 * physical area weighting and consistent nodal distribution.
 *
 * @param model_data Global nodal positions required for surface geometry.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void PLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // Surface coordinate mappings and normals require the global position field.
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "PLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    // Pressure is scalar, so its temporal multiplier can be evaluated once and
    // folded into the nominal magnitude before any geometric integration.
    const Precision scale           = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    const Precision scaled_pressure = pressure_ * scale;

    for (const ID surf_id : *region_) {
        // Resolve each selected surface and ignore absent entries safely.
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        // The surface performs quadrature and nodal shape-function weighting.
        // The callback supplies the pressure traction at the current global
        // integration-point position.
        surface->integrate_vector_field(
            node_positions,
            bc,
            [&](const Vec3& position) {
                // Recover the surface-local coordinates needed by the normal
                // implementation and apply pressure opposite to that normal.
                const Vec2 local = surface->global_to_local(position, node_positions);
                return -scaled_pressure * surface->normal(node_positions, local);
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the pressure load.
 *
 * The result identifies the target surface region, its size and the nominal
 * scalar pressure before temporal scaling.
 *
 * @return Human-readable load description.
 */
std::string PLoad::str() const {
    std::ostringstream os;

    // Report the surface set, its current size and the unscaled pressure.
    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", p=" << pressure_;

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}
} // namespace bc
} // namespace fem
