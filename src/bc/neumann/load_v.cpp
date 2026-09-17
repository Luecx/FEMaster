/**
 * @file load_v.cpp
 * @brief Implements consistent integration of distributed body-force densities.
 *
 * `VLoad` converts a prescribed force density to equivalent nodal forces through
 * the element weak-form contribution
 *
 *     f_e = integral_Omega_e N^T b dOmega.
 *
 * The structural element owns interpolation, quadrature and the physical volume
 * Jacobian. This implementation supplies the vector field `b(x)`, including
 * optional amplitude scaling and local-to-global transformation.
 *
 * The stored vector already has units of force per volume. Consequently the
 * element integration is called with density scaling disabled; multiplying by
 * material density here would incorrectly turn the load into a mass-dependent
 * inertia term.
 *
 * @see VLoad
 * @see Neumann
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#include "load_v.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {

/**
 * Integrates the prescribed body-force density over all selected elements.
 *
 * `NaN` components are replaced by zero before any arithmetic. The optional
 * amplitude is evaluated once because it scales the complete field uniformly.
 * For a global definition,
 *
 *     b(x) = a(t) b_0,
 *
 * while an oriented load uses
 *
 *     b(x) = a(t) A(x) b_local.
 *
 * The local basis is evaluated at each integration point so curvilinear
 * coordinate systems are represented correctly. The element integrator receives
 * `scale_by_density = false`, leaving only the geometric volume measure in the
 * finite-element integral.
 *
 * @param model_data Compiled element storage and nodal geometry.
 * @param rhs Generalized nodal RHS field modified in place.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Apply the nominal force density with unit amplitude
 *                         when true.
 */
void VLoad::apply(model::ModelData& model_data, model::Field& rhs, Precision time, bool ignore_amplitude) {
    // Validate the model geometry and semantic target required for volume
    // integration
    logging::error(model_data.positions != nullptr,
        "VLOAD: positions field is not initialized");
    logging::error(region_ != nullptr,
        "VLOAD: target element region is not initialized");

    // Convert the sparse input convention into a complete numerical vector
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

    auto [body_force_local, active] = sanitize_vector(values_);
    if (!active) {
        return;
    }

    // Apply the common scalar time history once to the nominal vector
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    body_force_local *= scale;

    // Integrate the body-force density independently over every selected
    // structural element
    for (ID element_id : *region_) {
        auto& element = model_data.elements[element_id];
        if (!element) {
            continue;
        }

        auto* structural = element->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        if (!orientation_) {
            // Global body-force density is constant over the element. The
            // structural integrator forms integral N^T b dOmega and scatters it
            // to the supplied global RHS field.
            auto body_force = [body_force_local](const Vec3&) -> Vec3 {
                return body_force_local;
            };

            structural->integrate_vector_field(rhs, false, body_force);
            continue;
        }

        // Evaluate A(x) at every quadrature point because a local coordinate
        // system may vary spatially over the element
        auto* orientation = orientation_.get();
        auto body_force = [orientation, body_force_local](const Vec3& position) -> Vec3 {
            const Vec3 local_point = orientation->to_local(position);
            const auto axes        = orientation->get_axes(local_point);
            return axes * body_force_local;
        };

        // VLOAD is a force density, not an acceleration. Density scaling must
        // therefore remain disabled in the element integration measure.
        structural->integrate_vector_field(rhs, false, body_force);
    }
}

/**
 * Builds a compact diagnostic representation of the body-force density.
 *
 * The output reports the semantic element target, nominal vector and optional
 * orientation and amplitude. The values are shown before coordinate
 * transformation and temporal scaling.
 *
 * @return Human-readable volume-load description.
 */
std::string VLoad::str() const {
    std::ostringstream os;

    os << "VLOAD: target=ELSET "
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
