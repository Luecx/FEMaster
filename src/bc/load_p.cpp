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

namespace {

inline Precision eval_scale(const Amplitude::Ptr& amp, Precision time) {
    return amp ? amp->evaluate(time) : Precision(1);
}

} // namespace

/**
 * @copydoc PLoad::apply
 */
void PLoad::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    logging::error(model_data.positions != nullptr, "positions field not set in model data");
    const auto& node_positions = *model_data.positions;

    const Precision scale = eval_scale(amplitude, time);
    const Precision scaled_pressure = pressure * scale;

    for (auto& surf_id : *region) {
        model_data.surfaces[surf_id]->apply_pload(node_positions, bc, scaled_pressure);
    }
}

/**
 * @copydoc PLoad::str
 */
std::string PLoad::str() const {
    std::ostringstream os;
    os << "PLOAD: target=SFSET " << (region ? region->name : std::string("?"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", p=" << pressure;
    if (amplitude) {
        os << ", amplitude=" << amplitude->name;
    }
    return os.str();
}

} // namespace bc
} // namespace fem
