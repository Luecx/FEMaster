/**
 * @file load_v.cpp
 * @brief Implements density-scaled distributed body-force assembly.
 *
 * The body-load vector is sanitized, optionally amplitude-scaled and supplied
 * as a spatial vector field to each selected structural element. The element
 * formulation owns the integration measure, shape functions and material
 * density scaling. A local orientation is evaluated at every integration point
 * when its basis may vary in space.
 *
 * @see load_v.h
 * @see model::StructuralElement
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_v.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {
namespace {

/**
 * Converts a partially specified body-load vector into numerical components.
 *
 * NaN marks omitted components in the physical load definition. Those entries
 * become zero before amplitude scaling or element integration, while the active
 * flag identifies a completely unspecified vector that can be skipped.
 *
 * @param vec Nominal vector with optional NaN markers.
 * @return Sanitized vector and flag indicating at least one active component.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    // Inspect all components independently because the input may prescribe only
    // selected global or local directions.
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            vec[i] = Precision(0);
        } else {
            active = true;
        }
    }

    return {vec, active};
}

} // namespace

/**
 * Assembles distributed body force over the selected structural elements.
 *
 * Sparse input is sanitized and optionally amplitude-scaled. Structural element
 * formulations perform density-aware integration and consistent nodal
 * distribution. When a local coordinate system is assigned, its basis is
 * evaluated at each integration point before the vector enters element
 * integration.
 *
 * @param model_data Global fields and element data required for assembly.
 * @param rhs Generalized nodal right-hand-side field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void VLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    // Validate the model geometry, selected element region and structural target
    // field before constructing the numerical body-load vector.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "VLoad: target element region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "VLoad: target field must be a NODE field with at least six components");

    // Convert sparse component input into a complete vector and skip a wholly
    // unspecified load before touching the element container.
    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) return;

    // Apply the optional temporal multiplier once to the nominal vector because
    // the same amplitude scales every selected element.
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    local_values *= scale;

    // Delegate physical volume, area or line integration to every selected
    // element that implements structural behavior.
    for (const ID element_id : *region_) {
        // Resolve the compiled element pointer and tolerate empty storage slots.
        auto& element = model_data.elements[element_id];
        if (!element) continue;

        // Non-structural elements do not provide the density-scaled vector-field
        // integration required by a structural body load.
        auto structural = element->as<model::StructuralElement>();
        if (!structural) continue;

        if (!orientation_) {
            // In the global case the vector field is constant over the complete
            // element. The element supplies density, interpolation and geometry.
            structural->integrate_vector_field(
                rhs,
                true,
                [local_values](const Vec3&) -> Vec3 { return local_values; }
            );
            continue;
        }

        // Preserve the shared orientation for the duration of synchronous
        // element integration and transform the nominal local vector at every
        // physical integration-point position.
        auto* orientation = orientation_.get();
        structural->integrate_vector_field(
            rhs,
            true,
            [orientation, local_values](const Vec3& position) -> Vec3 {
                const Vec3 local_point = orientation->to_local(position);
                return orientation->get_axes(local_point) * local_values;
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the distributed body load.
 *
 * The result identifies the target element region, its size and the nominal
 * density-scaled vector before amplitude scaling or coordinate transformation.
 *
 * @return Human-readable load description.
 */
std::string VLoad::str() const {
    std::ostringstream os;

    // Print the target element set, its size and the nominal unscaled vector.
    os << "VLOAD: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", values=[" << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";

    // Append optional modifiers only when they are assigned.
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
