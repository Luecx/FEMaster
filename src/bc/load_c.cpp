/**
 * @file load_c.cpp
 * @brief Implements the concentrated nodal load boundary condition.
 *
 * @see src/bc/load_c.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_c.h"

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

std::string vec6_to_string(const Vec6& v) {
    std::ostringstream os;
    os << v[0] << ", " << v[1] << ", " << v[2] << ", "
       << v[3] << ", " << v[4] << ", " << v[5];
    return os.str();
}

} // namespace

/**
 * @copydoc CLoad::apply
 */
void CLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    const auto& node_positions = *model_data.positions;

    const Precision scale = eval_scale(amplitude, time);

    for (auto& node_id : *region) {
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));

        auto [force_local, force_active]   = sanitize_vector(values.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values.tail<3>());

        force_local *= scale;
        moment_local *= scale;

        if (!orientation) {
            if (force_active) {
                for (int i = 0; i < 3; ++i) {
                    bc(node_id, i) += force_local[i];
                }
            }
            if (moment_active) {
                for (int i = 0; i < 3; ++i) {
                    bc(node_id, static_cast<Dim>(i + 3)) += moment_local[i];
                }
            }
            continue;
        }

        const Vec3 local_point = orientation->to_local(position);
        const auto axes = orientation->get_axes(local_point);

        if (force_active) {
            const Vec3 global_force = axes * force_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, i) += global_force[i];
            }
        }

        if (moment_active) {
            const Vec3 global_moment = axes * moment_local;
            for (int i = 0; i < 3; ++i) {
                bc(node_id, static_cast<Dim>(i + 3)) += global_moment[i];
            }
        }
    }
}

/**
 * @copydoc CLoad::str
 */
std::string CLoad::str() const {
    std::ostringstream os;
    os << "CLOAD: target=NSET " << (region ? region->name : std::string("?"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", values=[" << vec6_to_string(values) << "]";
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
