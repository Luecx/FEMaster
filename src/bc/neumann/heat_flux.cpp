/**
 * @file heat_flux.cpp
 * @brief Implements prescribed surface heat-flux boundary conditions.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "heat_flux.h"

#include "../../core/logging.h"
#include "../../model/model_data.h"

#include <cmath>
#include <sstream>

namespace fem::bc {

void HeatFlux::apply(model::ModelData& model_data,
                     model::Field&     bc,
                     Precision         time,
                     bool              ignore_amplitude) {
    logging::error(region_ != nullptr,
        "HEATFLUX: target surface region is not set");
    logging::error(model_data.positions_reference != nullptr,
        "HEATFLUX: reference positions are not initialized");
    logging::error(bc.domain == model::FieldDomain::NODE && bc.components >= 3,
        "HEATFLUX: target field must be a NODE field with at least three components");
    logging::error(std::isfinite(heat_flux_),
        "HEATFLUX: prescribed heat flux must be finite");

    const Precision scale = amplitude_ && !ignore_amplitude
        ? amplitude_->evaluate(time)
        : Precision(1);
    const Vec3 value{heat_flux_ * scale, Precision(0), Precision(0)};

    for (ID surface_id : *region_) {
        logging::error(surface_id >= 0
                    && static_cast<Index>(surface_id) < model_data.surfaces.size(),
            "HEATFLUX: surface ", surface_id, " is outside the compiled surface domain");

        const auto& surface = model_data.surfaces[static_cast<Index>(surface_id)];
        logging::error(surface != nullptr,
            "HEATFLUX: surface ", surface_id, " is not initialized");

        surface->integrate_vector_field(
            *model_data.positions_reference,
            bc,
            [value](const Vec3&) -> Vec3 { return value; }
        );
    }
}

std::string HeatFlux::str() const {
    std::ostringstream os;
    os << "HEATFLUX: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", value=" << heat_flux_;

    if (amplitude_)
        os << ", amplitude=" << amplitude_->name;

    return os.str();
}

} // namespace fem::bc
