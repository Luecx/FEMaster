/**
 * @file load_p.cpp
 * @brief Implements the uniform surface pressure load.
 *
 * @see src/bc/load_p.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "load_p.h"

#include "../core/logging.h"
#include "../model/model_data.h"

#include <sstream>

namespace fem {
namespace bc {
/**
 * @copydoc PLoad::apply
 */
void PLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    logging::error(region_ != nullptr, "PLoad: target surface region not set");
    const auto& node_positions = *model_data.positions;

    // Pressure is scalar, so amplitude scaling can be applied once before the
    // geometry-specific surface integration distributes it to nodes.
    const Precision scale           = amplitude_ ? amplitude_->evaluate(time) : Precision(1);
    const Precision scaled_pressure = pressure_ * scale;

    for (const ID surf_id : *region_) {
        model_data.surfaces[surf_id]->apply_pload(node_positions, bc, scaled_pressure);
    }
}

/**
 * @copydoc PLoad::str
 */
std::string PLoad::str() const {
    std::ostringstream os;

    os << "PLOAD: target=SFSET "
       << (region_ ? region_->name : std::string("?"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", p=" << pressure_;

    if (amplitude_) {
        os << ", amplitude=" << amplitude_->name;
    }

    return os.str();
}
} // namespace bc
} // namespace fem
