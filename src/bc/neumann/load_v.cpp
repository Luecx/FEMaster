/**
 * @file load_v.cpp
 * @brief Implements density-scaled distributed body-force assembly.
 *
 * The body-load vector is sanitized, optionally amplitude-scaled and supplied
 * as a spatial vector field to each selected structural element. The element
 * formulation owns the integration measure, shape functions and material
 * density scaling.
 *
 * @see load_v.h
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
 * Converts a partially specified vector into an assembly-safe vector.
 *
 * @param vec Vector containing finite values and optional NaN markers.
 * @return Sanitized vector and an active-component flag.
 */
std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;

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
 * @param model_data Global fields and element data required for assembly.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void VLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // Validate the model geometry and selected element region.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "VLoad: target element region not set");

    // Sanitize sparse input and apply the common amplitude once.
    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    local_values *= scale;

    for (const ID el_id : *region_) {
        auto& element = model_data.elements[el_id];
        if (!element) {
            continue;
        }

        auto structural = element->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        if (!orientation_) {
            const Vec3 value = local_values;
            structural->integrate_vector_field(
                bc,
                true,
                [value](const Vec3&) -> Vec3 { return value; }
            );
            continue;
        }

        // Transform the local vector independently at each integration point.
        auto* orientation = orientation_.get();
        structural->integrate_vector_field(
            bc,
            true,
            [orientation, local_values](const Vec3& x) -> Vec3 {
                const Vec3 local_point = orientation->to_local(x);
                const auto axes        = orientation->get_axes(local_point);
                return axes * local_values;
            }
        );
    }
}

/**
 * Builds the diagnostic representation of the body load.
 *
 * @return Human-readable load description.
 */
std::string VLoad::str() const {
    std::ostringstream os;
    os << "VLOAD: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", values=[" << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";

    if (orientation_)
        os << ", orientation=" << orientation_->name;
    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
