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
/**
 * @brief Converts a partially specified traction vector into usable values.
 *
 * `NaN` means an unset component in user input. The assembly routines expect a
 * complete vector, so unset components are converted to zero while the returned
 * flag records whether the load has any active component at all.
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
 * @copydoc DLoad::apply
 */
void DLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "DLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude_ ? amplitude_->evaluate(time) : Precision(1);
    local_values *= scale;

    for (const ID surf_id : *region_) {
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        if (!orientation_) {
            surface->apply_dload(node_positions, bc, local_values);
            continue;
        }

        // For oriented tractions the local basis can vary over the surface.
        // Therefore every nodal contribution is rotated at its own node
        // position before it is added to the global RHS.
        const auto contributions = surface->shape_function_integral(node_positions);
        int        local_idx     = 0;
        for (auto node_it = surface->begin(); node_it != surface->end(); ++node_it, ++local_idx) {
            const ID node_id = *node_it;

            const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));
            const Vec3 local_point = orientation_->to_local(position);
            const auto axes = orientation_->get_axes(local_point);
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

    os << "DLOAD: target=SFSET "
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
