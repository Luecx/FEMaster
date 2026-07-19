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
void DLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "DLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) {
        return;
    }

    const Precision scale = amplitude_ && !ignore_amplitude ? amplitude_->evaluate(time) : Precision(1);
    local_values *= scale;

    for (const ID surf_id : *region_) {
        auto surface = model_data.surfaces[surf_id];
        if (!surface) {
            continue;
        }

        if (!orientation_) {
            surface->integrate_vector_field(
                node_positions,
                bc,
                [&](const Vec3&) -> Vec3 {
                    return local_values;
                }
            );
            continue;
        }

        // Rotate oriented tractions at each integration point
        surface->integrate_vector_field(
            node_positions,
            bc,
            [&](const Vec3& position) -> Vec3 {
                const Vec3 local_point = orientation_->to_local(position);
                const auto axes        = orientation_->get_axes(local_point);
                return axes * local_values;
            }
        );
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
