/**
 * @file load_collector.cpp
 * @brief Implements load collector construction and RHS dispatch.
 */

#include "load_collector.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

void LoadCollector::apply(model::ModelData& model_data,
                          model::Field&     rhs,
                          Precision         time) {
    for (const auto& load : this->_data) {
        if (load) load->apply(model_data, rhs, time, false);
    }
}

} // namespace fem::bc
