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
/**
 * @brief Converts a partially specified 3-vector into an assemble-ready vector.
 *
 * A `NaN` component means "not specified" in the load input. Assembly code is
 * simpler and safer when those entries are replaced by zero and the helper also
 * reports whether at least one component was actually prescribed.
 *
 * @param vector Input vector that may contain `NaN` markers.
 * @return std::pair<Vec3, bool> Sanitized vector and an active-component flag.
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
 * @copydoc CLoad::apply
 */
void CLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "CLoad: target node region not set");
    const auto& node_positions = *model_data.positions;

    const Precision scale = amplitude_ ? amplitude_->evaluate(time) : Precision(1);

    for (const ID node_id : *region_) {
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));

        auto [force_local , force_active]  = sanitize_vector(values_.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());

        // The amplitude is applied after NaN sanitizing so unspecified
        // components stay exactly zero for both force and moment parts.
        force_local *= scale;
        moment_local *= scale;

        if (!orientation_) {
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

        // Local loads are interpreted in the orientation basis at the current
        // node position. Translational and rotational components use the same
        // three axes, but are written to different RHS columns.
        const Vec3 local_point = orientation_->to_local(position);
        const auto axes        = orientation_->get_axes(local_point);

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

    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", values=["
       << values_[0] << ", "
       << values_[1] << ", "
       << values_[2] << ", "
       << values_[3] << ", "
       << values_[4] << ", "
       << values_[5]
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
