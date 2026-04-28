/**
 * @file load_v.cpp
 * @brief Implements the volumetric element load.
 *
 * @see src/bc/load_v.h
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
 * @brief Converts a partially specified body-force vector into usable values.
 *
 * `NaN` denotes an unset component. The element integration callback expects a
 * complete vector, so unset entries are converted to zero and the returned flag
 * allows callers to skip a fully empty load.
 *
 * @param vector Input vector that may contain `NaN` markers.
 * @return std::pair<Vec3, bool> Sanitized vector and active-component flag.
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
 * @copydoc VLoad::apply
 */
void VLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "VLoad: target element region not set");

    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude_ ? amplitude_->evaluate(time) : Precision(1);
    local_values *= scale;

    for (const ID el_id : *region_) {
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        if (!orientation_) {
            const Vec3 f0 = local_values;
            // Without an orientation the body force is constant in global
            // coordinates and can be integrated directly over the element.
            auto f = [f0](const Vec3& /*x*/) -> Vec3 { return f0; };
            structural->integrate_vec_field(bc, /*scale_by_density=*/true, f);
        } else {
            const Vec3 f_local = local_values;
            auto*      ori     = orientation_.get();
            // Oriented body loads may use a spatially varying local basis. The
            // callback rotates the local vector at each integration point.
            auto f = [ori, f_local](const Vec3& x) -> Vec3 {
                const Vec3 local_point = ori->to_local(x);
                const auto axes        = ori->get_axes(local_point);
                return axes * f_local;
            };
            structural->integrate_vec_field(bc, /*scale_by_density=*/true, f);
        }
    }
}

/**
 * @copydoc VLoad::str
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
} // namespace bc
} // namespace fem
