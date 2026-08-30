/**
 * @file load_d.cpp
 * @brief Implements consistent distributed surface traction assembly.
 */

#include "load_d.h"

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

void DLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "DLoad: target surface region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "DLoad: target field must be a NODE field with at least six components");

    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) return;

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    local_values *= scale;

    const auto& positions = *model_data.positions;
    for (const ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        if (!orientation_) {
            surface->integrate_vector_field(
                positions,
                rhs,
                [local_values](const Vec3&) -> Vec3 { return local_values; }
            );
            continue;
        }

        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                const Vec3 local_point = orientation_->to_local(position);
                return orientation_->get_axes(local_point) * local_values;
            }
        );
    }
}

std::string DLoad::str() const {
    std::ostringstream os;
    os << "DLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", values=[" << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;
    return os.str();
}

} // namespace fem::bc
