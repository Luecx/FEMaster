/**
 * @file support_collector.cpp
 * @brief Implements collective support equation generation.
 *
 * @author Finn Eggers
 * @date 19.08.2026
 */

#include "support_collector.h"

namespace fem::bc {

SupportCollector::SupportCollector(const std::string& name)
    : model::Collection<Support>(name, true, false) {}

constraint::Equations SupportCollector::get_equations(model::ModelData& model_data) {
    constraint::Equations equations{};
    for (Support& support : this->_data) {
        support.apply(model_data, equations);
    }
    return equations;
}

} // namespace fem::bc
