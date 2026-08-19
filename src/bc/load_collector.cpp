/**
 * @file load_collector.cpp
 * @brief Implements collective application of load objects.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "load_collector.h"

#include "../data/field.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Load::Ptr>(name) {}

void LoadCollector::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    for (const auto& load : this->_data) {
        if (load) {
            load->apply(model_data, bc, time);
        }
    }
}

} // namespace fem::bc
