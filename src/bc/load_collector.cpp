/**
 * @file load_collector.cpp
 * @brief Implements structural load collector dispatch.
 */

#include "load_collector.h"

#include "neumann/neumann.h"
#include "../core/logging.h"

#include <memory>

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

void LoadCollector::apply(model::ModelData& model_data,
                          model::Field&     rhs,
                          Precision         time) {
    for (const auto& load : this->_data) {
        if (!load) continue;

        if (auto neumann = std::dynamic_pointer_cast<StructuralNeumann>(load)) {
            neumann->apply(model_data, rhs, time, false);
            continue;
        }

        logging::error(false,
            "Load collector ", name,
            " contains a non-structural condition ", load->str());
    }
}

} // namespace fem::bc
