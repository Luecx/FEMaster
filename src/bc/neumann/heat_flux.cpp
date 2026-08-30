/**
 * @file heat_flux.cpp
 * @brief Implements consistent surface heat-flux assembly.
 */

#include "heat_flux.h"
#include "thermal_surface.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <sstream>

namespace fem::bc {

void HeatFlux::apply(model::ModelData&   model_data,
                     const SystemDofIds& thermal_dof_ids,
                     TripletList&        matrix_terms,
                     DynamicVector&      heat_source) const {
    (void) matrix_terms;

    logging::error(region_ != nullptr,
        "HeatFlux: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HeatFlux: reference positions field is not initialized");

    const auto& positions = *model_data.positions_reference;

    for (const ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < static_cast<Index>(model_data.surfaces.size()),
            "HeatFlux: surface id ", surface_id, " is out of range");

        const auto& surface = model_data.surfaces[static_cast<std::size_t>(surface_id)];
        logging::error(surface != nullptr,
            "HeatFlux: surface ", surface_id, " is not initialized");

        detail::visit_thermal_surface(*surface, [&](const auto& typed_surface) {
            detail::assemble_flux_on_surface(
                typed_surface,
                positions,
                thermal_dof_ids,
                flux_,
                heat_source
            );
        });
    }
}

std::string HeatFlux::str() const {
    std::ostringstream os;
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", q=" << flux_;
    return os.str();
}

} // namespace fem::bc
