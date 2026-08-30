/**
 * @file load_p.cpp
 * @brief Implements pressure integration along surface normals.
 */

#include "load_p.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

void PLoad::apply(model::ModelData& model_data,
                  model::Field&     rhs,
                  Precision         time,
                  bool              ignore_amplitude) {
    logging::error(model_data.positions != nullptr,
        "positions field not set in model data");
    logging::error(region_ != nullptr,
        "PLoad: target surface region not set");
    logging::error(rhs.domain == model::FieldDomain::NODE && rhs.components >= 6,
        "PLoad: target field must be a NODE field with at least six components");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Precision pressure = pressure_ * scale;
    const auto& positions = *model_data.positions;

    for (const ID surface_id : *region_) {
        const auto& surface = model_data.surfaces[surface_id];
        if (!surface) continue;

        surface->integrate_vector_field(
            positions,
            rhs,
            [&](const Vec3& position) -> Vec3 {
                const Vec2 local = surface->global_to_local(position, positions);
                return -pressure * surface->normal(positions, local);
            }
        );
    }
}

std::string PLoad::str() const {
    std::ostringstream os;
    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", p=" << pressure_;
    if (amplitude_) os << ", amplitude=" << amplitude_->name;
    return os.str();
}

} // namespace fem::bc
