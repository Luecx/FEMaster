/**
 * @file load_c.cpp
 * @brief Implements concentrated nodal force and moment assembly.
 */

#include "load_c.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>
#include <utility>

namespace fem::bc {
namespace {

std::pair<Vec3, bool> sanitize_vector(Vec3 vec) {
    bool active = false;
    for (int i = 0; i < 3; ++i) {
        if (std::isnan(vec[i])) vec[i] = Precision(0);
        else active = true;
    }
    return {vec, active};
}

} // namespace

void CLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "CLoad: target node region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "CLoad: target field must be a NODE field with at least six components");

    const auto& node_positions = *model_data.positions;
    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);

    for (const ID node_id : *region_) {
        const Vec3 position = node_positions.row_vec3(static_cast<Index>(node_id));
        auto [force_local, force_active] = sanitize_vector(values_.head<3>());
        auto [moment_local, moment_active] = sanitize_vector(values_.tail<3>());
        force_local  *= scale;
        moment_local *= scale;

        if (!orientation_) {
            if (force_active) {
                for (Dim i = 0; i < 3; ++i) rhs(node_id, i) += force_local[i];
            }
            if (moment_active) {
                for (Dim i = 0; i < 3; ++i) rhs(node_id, i + 3) += moment_local[i];
            }
            continue;
        }

        const Vec3 local_point = orientation_->to_local(position);
        const auto axes = orientation_->get_axes(local_point);

        if (force_active) {
            const Vec3 global_force = axes * force_local;
            for (Dim i = 0; i < 3; ++i) rhs(node_id, i) += global_force[i];
        }
        if (moment_active) {
            const Vec3 global_moment = axes * moment_local;
            for (Dim i = 0; i < 3; ++i) rhs(node_id, i + 3) += global_moment[i];
        }
    }
}

std::string CLoad::str() const {
    std::ostringstream os;
    os << "CLOAD: target=NSET "
       << (region_ ? region_->name : std::string("?")) << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0) << "), values=["
       << values_[0] << ", " << values_[1] << ", " << values_[2] << ", "
       << values_[3] << ", " << values_[4] << ", " << values_[5] << "]";
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;
    return os.str();
}

} // namespace fem::bc
