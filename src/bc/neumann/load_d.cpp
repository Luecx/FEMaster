/**
 * @file load_d.cpp
 * @brief Implements consistent integration of distributed surface tractions.
 *
 * For every selected surface, the continuous traction field is converted to
 * equivalent nodal forces through the standard weak-form contribution
 *
 *     f_e = integral_Gamma_e N^T t dGamma.
 *
 * The surface implementation owns interpolation, quadrature and the physical
 * area Jacobian. `DLoad` provides the traction integrand, including optional
 * amplitude scaling and local-to-global coordinate transformation. Curvilinear
 * coordinate systems are evaluated independently at every integration point.
 *
 * @see DLoad
 * @see Neumann
 * @see model::SurfaceInterface
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_d.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {

/**
 * Integrates the prescribed traction over all selected surfaces.
 *
 * `NaN` input components are first replaced by zero. The common temporal
 * amplitude is then applied to the nominal vector. For a global traction the
 * integrand is spatially constant,
 *
 *     t(x) = a(t) t_0,
 *
 * whereas an oriented traction uses
 *
 *     t(x) = a(t) A(x) t_local.
 *
 * `SurfaceInterface::integrate_vector_field()` supplies the interpolation and
 * physical surface measure, so this routine only defines the vector field being
 * integrated and selects the target surfaces.
 *
 * @param model_data Global nodal geometry and compiled surface storage.
 * @param rhs Generalized nodal RHS field modified in place.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal traction with unit amplitude when true.
 */
void DLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // Validate the geometric data and semantic target required for surface
    // integration
    logging::error(model_data.positions != nullptr,
        "DLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "DLOAD: target surface region is not initialized");

    const auto& node_positions = *model_data.positions;

    // Convert the sparse component convention into a numerical traction vector
    auto sanitize_vector = [](Vec3 vector) {
        bool active = false;
        for (Dim component = 0; component < 3; ++component) {
            if (std::isnan(vector[component])) {
                vector[component] = Precision(0);
            } else {
                active = true;
            }
        }
        return std::pair<Vec3, bool>{vector, active};
    };

    auto [traction_local, active] = sanitize_vector(values_);
    if (!active) {
        return;
    }

    // Apply the scalar time history once because it is common to every surface
    // and every quadrature point
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    traction_local *= scale;

    // Integrate the traction contribution independently on every selected surface
    for (ID surface_id : *region_) {
        auto surface = model_data.surfaces[surface_id];
        if (!surface) {
            continue;
        }

        if (!orientation_) {
            // The traction is constant in the global basis. The surface routine
            // evaluates N^T t and multiplies by the physical quadrature measure.
            surface->integrate_vector_field(
                node_positions,
                rhs,
                [&](const Vec3&) -> Vec3 {
                    return traction_local;
                }
            );
            continue;
        }

        // A local traction must be transformed inside the integrand because the
        // basis A(x) may vary over a curved or spatially defined coordinate system.
        surface->integrate_vector_field(
            node_positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                const Vec3 local_point = orientation_->to_local(position);
                const auto axes        = orientation_->get_axes(local_point);
                return axes * traction_local;
            }
        );
    }
}

/**
 * Builds a compact diagnostic representation of the distributed traction.
 *
 * The output reports the semantic surface target, the nominal traction vector
 * before amplitude scaling and coordinate transformation, and the optional
 * modifier names.
 *
 * @return Human-readable traction-load description.
 */
std::string DLoad::str() const {
    std::ostringstream os;

    os << "DLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", values=["
       << values_[0] << ", "
       << values_[1] << ", "
       << values_[2]
       << "]";

    if (orientation_) {
        os << ", orientation=" << orientation_->name;
    }

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}

} // namespace fem::bc
