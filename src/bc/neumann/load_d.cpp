/**
 * @file load_d.cpp
 * @brief Implements consistent integration of distributed surface tractions.
 *
 * The implementation sanitizes partially specified traction components,
 * applies optional temporal scaling and delegates geometric quadrature and
 * nodal distribution to each selected surface. Local traction vectors are
 * rotated into the global basis independently at every integration point.
 *
 * @see load_d.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_d.h"

#include "../../core/logging.h"
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
 * Integrates distributed traction over the selected surfaces.
 *
 * @param model_data Global fields required for surface integration.
 * @param bc Generalized nodal field receiving the contribution.
 * @param time Analysis time used for amplitude evaluation.
 * @param ignore_amplitude Whether amplitude scaling is disabled.
 */
void DLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time, bool ignore_amplitude) {
    // Validate surface geometry and the target region.
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "DLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    // Sanitize sparse input and apply the common amplitude once.
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
                [&](const Vec3&) -> Vec3 { return local_values; }
            );
            continue;
        }

        // Evaluate position-dependent local axes at every surface quadrature point.
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
 * Builds the diagnostic representation of the distributed traction load.
 *
 * @return Human-readable load description.
 */
std::string DLoad::str() const {
    std::ostringstream os;
    os << "DLOAD: target=SFSET "
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
