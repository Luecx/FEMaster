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

inline Precision eval_scale(const Amplitude::Ptr& amp, Precision time) {
    return amp ? amp->evaluate(time) : Precision(1);
}

std::string vec3_to_string(const Vec3& v) {
    std::ostringstream os;
    os << v[0] << ", " << v[1] << ", " << v[2];
    return os.str();
}

} // namespace

/**
 * @copydoc VLoad::apply
 */
void VLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");

    auto [local_values, has_values] = sanitize_vector(values);
    if (!has_values) {
        return;
    }

    const Precision scale = eval_scale(amplitude, time);
    local_values *= scale;

    for (auto& el_id : *region) {
        auto& el_ptr = model_data.elements[el_id];
        if (!el_ptr) {
            continue;
        }

        auto structural = el_ptr->as<model::StructuralElement>();
        if (!structural) {
            continue;
        }

        if (!orientation) {
            const Vec3 f0 = local_values;
            auto f = [f0](const Vec3& /*x*/) -> Vec3 { return f0; };
            structural->integrate_vec_field(bc, /*scale_by_density=*/true, f);
        } else {
            const Vec3 f_local = local_values;
            auto* ori = orientation.get();
            auto f = [ori, f_local](const Vec3& x) -> Vec3 {
                const Vec3 local_point = ori->to_local(x);
                const auto axes = ori->get_axes(local_point);
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
    os << "VLOAD: target=ELSET " << (region ? region->name : std::string("?"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", values=[" << vec3_to_string(values) << "]";
    if (orientation) {
        os << ", orientation=" << orientation->name;
    }
    if (amplitude) {
        os << ", amplitude=" << amplitude->name;
    }
    return os.str();
}

} // namespace bc
} // namespace fem
