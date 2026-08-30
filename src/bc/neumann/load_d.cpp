/**
 * @file load_d.cpp
 * @brief Implements consistent integration of distributed surface tractions.
 *
 * The implementation sanitizes partially specified traction components,
 * applies optional temporal scaling and delegates geometric quadrature and
 * nodal distribution to each selected surface. Local traction vectors are
 * rotated into the global basis independently at every integration point.
 *
 * @see load_d.h
 * @see SurfaceInterface
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_d.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {
namespace {

/**
 * Converts a partially specified traction vector into numerical components.
 *
 * NaN marks an omitted component in the load definition. Omitted components are
 * replaced by zero before amplitude scaling and surface integration, while the
 * returned flag records whether the traction contains any prescribed component.
 *
 * @param vec Nominal traction vector with optional NaN markers.
 * @return Sanitized vector and flag indicating at least one active component.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    // Inspect each component independently because a traction may be prescribed
    // only in selected directions.
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            // Omitted components must not contaminate integration arithmetic.
            vec[i] = Precision(0);
        } else {
            active = true;
        }
    }

    return {vec, active};
}

} // namespace

/**
 * Integrates distributed traction over the selected surfaces.
 *
 * The traction is sanitized and optionally amplitude-scaled. Surface geometry
 * supplies quadrature, physical area weighting and consistent nodal
 * distribution. When an orientation is assigned, the nominal local traction is
 * transformed into global Cartesian coordinates separately at every quadrature
 * point so position-dependent bases remain valid.
 *
 * @param model_data Global fields required for surface integration.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void DLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    // Validate surface geometry, the target region and the structural RHS layout
    // before constructing the numerical traction vector.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "DLoad: target surface region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "DLoad: target field must be a NODE field with at least six components");

    // Convert sparse component input into a complete numerical vector. A load
    // with no finite component can return before traversing any surfaces.
    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) return;

    // The amplitude is common to the complete surface region and is therefore
    // evaluated only once per application.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    local_values *= scale;

    // Integrate the traction independently over every selected compiled surface.
    const auto& positions = *model_data.positions;
    for (const ID surface_id : *region_) {
        // Resolve the surface by its global identifier. Empty compiled slots do
        // not contribute to the load field.
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        if (!orientation_) {
            // A global traction is spatially constant. The surface integrator
            // evaluates the callback at quadrature points, multiplies by shape
            // functions and the physical Jacobian, and scatters nodal forces.
            surface->integrate_vector_field(
                positions,
                rhs,
                [local_values](const Vec3&) -> Vec3 { return local_values; }
            );
            continue;
        }

        // A local traction must be transformed where it is integrated because
        // coordinate-system axes may vary over the physical surface.
        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                // Evaluate the basis in the coordinate system's local
                // parameterization and map the nominal local vector to global
                // Cartesian components.
                const Vec3 local_point = orientation_->to_local(position);
                return orientation_->get_axes(local_point) * local_values;
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the distributed traction load.
 *
 * The result identifies the target surface region, its size and the nominal
 * traction components before amplitude scaling and coordinate transformation.
 *
 * @return Human-readable load description.
 */
std::string DLoad::str() const {
    std::ostringstream os;

    // Report the surface set, its size and all nominal traction components.
    os << "DLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", values=[" << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";

    // Append optional modifiers only when they are assigned.
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
