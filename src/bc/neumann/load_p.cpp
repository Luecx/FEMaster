/**
 * @file load_p.cpp
 * @brief Implements consistent pressure integration over finite-element surfaces.
 *
 * The scalar pressure is converted to the geometric traction
 *
 *     t(x) = -p n(x),
 *
 * where `n(x)` is the unit normal computed from the nodal geometry supplied by
 * `ModelData::positions`. Each selected surface then evaluates the consistent
 * weak-form contribution
 *
 *     f_e = -integral_Gamma_e N^T p n dGamma.
 *
 * Surface interpolation, quadrature and the physical surface Jacobian remain
 * encapsulated by the surface implementation. `PLoad` supplies only the scalar
 * magnitude, amplitude scaling and normal traction field.
 *
 * @see PLoad
 * @see Neumann
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_p.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

/**
 * Integrates the pressure traction over every selected surface.
 *
 * The optional amplitude is evaluated once and folded into the scalar pressure.
 * At each integration point, the surface position is mapped back to natural
 * coordinates so the geometric normal can be evaluated from the same nodal
 * configuration used by the surface integrator. The resulting traction
 *
 *     t = -a(t) p n
 *
 * is then weighted by the surface interpolation and physical area measure by
 * `integrate_vector_field()`.
 *
 * The normal follows the geometry represented by `model_data.positions`; this
 * routine does not independently choose between reference and current
 * configurations.
 *
 * @param model_data Nodal geometry and compiled surface storage.
 * @param rhs Generalized nodal RHS field modified in place.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal pressure with unit amplitude when true.
 */
void PLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // Validate the physical definition and geometry required for surface
    // integration
    logging::error(model_data.positions != nullptr,
        "PLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "PLOAD: target surface region is not initialized");
    logging::error(std::isfinite(pressure_),
        "PLOAD: pressure must be finite");

    const auto& node_positions = *model_data.positions;

    // Apply the common scalar time history before geometric integration
    const Precision scale           = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    const Precision scaled_pressure = pressure_ * scale;

    // Integrate the pressure traction independently on every selected surface
    for (ID surface_id : *region_) {
        auto surface = model_data.surfaces[surface_id];
        if (!surface) {
            continue;
        }

        // The surface routine constructs the consistent N^T t contribution and
        // physical quadrature measure. This callback provides only t(x).
        surface->integrate_vector_field(
            node_positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                // Recover the natural coordinates corresponding to this global
                // integration point and evaluate the unit normal from the same
                // nodal geometry used for integration.
                const Vec2 local  = surface->global_to_local(position, node_positions);
                const Vec3 normal = surface->normal(node_positions, local);

                return -scaled_pressure * normal;
            }
        );
    }
}

/**
 * Builds a compact diagnostic representation of the pressure definition.
 *
 * The output reports the semantic surface target, nominal scalar pressure and
 * optional amplitude. It does not evaluate surface normals or temporal scaling.
 *
 * @return Human-readable pressure-load description.
 */
std::string PLoad::str() const {
    std::ostringstream os;

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

} // namespace fem::bc
