/**
 * @file convection.cpp
 * @brief Implements film-convection matrix and source assembly.
 */

#include "convection.h"
#include "thermal_surface.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

void Convection::apply(model::ModelData&   model_data,
                       const SystemDofIds& thermal_dof_ids,
                       TripletList&        matrix_terms,
                       DynamicVector&      heat_source) const {
    logging::error(region_ != nullptr,
        "Convection: target surface region is not set");
    logging::error(film_coefficient_ >= Precision(0),
        "Convection: film coefficient must be non-negative");
    logging::error(model_data.positions_reference != nullptr,
        "Convection: reference positions field is not initialized");

    const auto& positions = *model_data.positions_reference;

    for (const ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
            "Convection: surface id ", surface_id, " is out of range");

        const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
        logging::error(surface != nullptr,
            "Convection: surface ", surface_id, " is not initialized");

        detail::visit_thermal_surface(*surface, [&](const auto& typed_surface) {
            detail::assemble_convection_on_surface(
                typed_surface,
                positions,
                thermal_dof_ids,
                film_coefficient_,
                ambient_temperature_,
                matrix_terms,
                heat_source
            );
        });
    }
}

std::string Convection::str() const {
    std::ostringstream os;
    os << "CONVECTION: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", h=" << film_coefficient_
       << ", T_inf=" << ambient_temperature_;
    return os.str();
}

} // namespace fem::bc
