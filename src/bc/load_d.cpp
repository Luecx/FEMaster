/**
 * @file load_d.cpp
 * @brief Implements the distributed surface traction load.
 *
 * @see src/bc/load_d.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_d.h"

#include "../core/logging.h"
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
 * @copydoc DLoad::apply
 */
void DLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    const auto& node_positions = *model_data.positions;

    auto [local_values, has_values] = sanitize_vector(values);
    if (!has_values) {
        return;
    }

    const Precision scale = eval_scale(amplitude, time);
    local_values *= scale;

    for (auto& surf_id : *region) {
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        if (!orientation) {
            surface->apply_dload(node_positions, bc, local_values);
            continue;
        }

        const auto contributions = surface->shape_function_integral(node_positions);
        int local_idx = 0;
        for (auto node_it = surface->begin(); node_it != surface->end(); ++node_it, ++local_idx) {
            const ID node_id = *node_it;

            const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));
            const Vec3 local_point = orientation->to_local(position);
            const auto axes = orientation->get_axes(local_point);
            const Vec3 global_values = axes * local_values;

            bc(node_id, 0) += contributions(local_idx) * global_values[0];
            bc(node_id, 1) += contributions(local_idx) * global_values[1];
            bc(node_id, 2) += contributions(local_idx) * global_values[2];
        }
    }
}

/**
 * @copydoc DLoad::str
 */
std::string DLoad::str() const {
    std::ostringstream os;
    os << "DLOAD: target=SFSET " << (region ? region->name : std::string("?"))
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
