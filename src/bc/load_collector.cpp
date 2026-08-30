/**
 * @file load_collector.cpp
 * @brief Implements collective application of Neumann conditions.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "load_collector.h"

#include "../data/field.h"

namespace fem::bc {

LoadCollector::LoadCollector(const std::string& name)
    : model::Collection<Neumann::Ptr>(name) {}

void LoadCollector::apply(model::ModelData& model_data, model::Field& bc, Precision time) {
    for (const auto& load : this->_data) {
        if (load) {
            load->apply(model_data, bc, time);
        }
    }
}

} // namespace fem::bc
