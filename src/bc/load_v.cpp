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
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_v.h"

#include "../core/logging.h"
#include "../model/element/element_structural.h"
#include "../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem {
namespace bc {
namespace {

/**
 * Converts a partially specified vector into an assembly-safe vector.
 *
 * Input decks use NaN to mark omitted components. Those entries become zero,
 * while the returned flag records whether any component was prescribed.
 *
 * @param vec Vector containing finite values and optional NaN markers.
 * @return The sanitized vector and an active-component flag.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) {
            vec[i] = 0.0;
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
 * Sparse input is sanitized and optionally amplitude-scaled. Element
 * formulations perform density-aware volume integration and nodal distribution;
 * an optional local basis is evaluated at each integration point.
 *
 * @param model_data Global fields and element data required for assembly.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void VLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // Position data is required by element integration and by local coordinate
    // systems evaluated at integration points.
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "VLoad: target element region not set");

    // Convert sparse component input into a complete vector and skip a wholly
    // unspecified load before touching the element container.
    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) {
        return;
    }

    // Apply the optional temporal multiplier once to the nominal vector.
    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    local_values *= scale;

    for (const ID el_id : *region_) {
        // Resolve the element pointer and tolerate empty slots in model storage.
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        // Only structural elements provide the vector-field integration needed
        // to convert a body load into equivalent nodal forces.
        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        if (!orientation_) {
            // In the global case the vector field is constant over the complete
            // element. Capture a value copy so the callback remains independent
            // of subsequent local variables.
            const Vec3 f0 = local_values;
            auto f = [f0](const Vec3& /*x*/) -> Vec3 { return f0; };

            // The element multiplies the field by material density, shape
            // functions and its geometric integration measure.
            structural->integrate_vector_field(bc, /*scale_by_density=*/true, f);
        } else {
            // Preserve the local nominal vector and a raw pointer to the shared
            // orientation for the duration of this synchronous integration call.
            const Vec3 f_local = local_values;
            auto*      ori     = orientation_.get();

            // Transform the local vector independently at each integration
            // point so spatially varying coordinate-system axes are respected.
            auto f = [ori, f_local](const Vec3& x) -> Vec3 {
                const Vec3 local_point = ori->to_local(x);
                const auto axes        = ori->get_axes(local_point);
                return axes * f_local;
            };
            structural->integrate_vector_field(bc, /*scale_by_density=*/true, f);
        }
    }
}

/**
 * Builds the diagnostic representation of the distributed body load.
 *
 * The result identifies the target element region, its size, nominal vector and
 * density-scaling configuration.
 *
 * @return Human-readable load description.
 */
std::string VLoad::str() const {
    std::ostringstream os;

    // Print the target element set, its size and the nominal unscaled vector.
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
} // namespace bc
} // namespace fem
