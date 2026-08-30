/**
 * @file load_v.cpp
 * @brief Implements density-scaled distributed body-force assembly.
 */

#include "load_v.h"

#include "../../core/logging.h"
#include "../../model/element/element_structural.h"
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

void VLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "VLoad: target element region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "VLoad: target field must be a NODE field with at least six components");

    auto [local_values, has_values] = sanitize_vector(values_);
    if (!has_values) return;

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    local_values *= scale;

    for (const ID element_id : *region_) {
        auto& element = model_data.elements[element_id];
        if (!element) continue;

        auto structural = element->as<model::StructuralElement>();
        if (!structural) continue;

        if (!orientation_) {
            structural->integrate_vector_field(
                rhs,
                true,
                [local_values](const Vec3&) -> Vec3 { return local_values; }
            );
            continue;
        }

        auto* orientation = orientation_.get();
        structural->integrate_vector_field(
            rhs,
            true,
            [orientation, local_values](const Vec3& position) -> Vec3 {
                const Vec3 local_point = orientation->to_local(position);
                return orientation->get_axes(local_point) * local_values;
            }
        );
    }
}

std::string VLoad::str() const {
    std::ostringstream os;
    os << "VLOAD: target=ELSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", values=[" << values_[0] << ", " << values_[1] << ", " << values_[2] << "]";
    if (orientation_) os << ", orientation=" << orientation_->name;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;
    return os.str();
}

} // namespace fem::bc
