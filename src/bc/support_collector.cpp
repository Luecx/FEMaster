/**
 * @file support_collector.cpp
 * @brief Implements construction and structural dispatch of support collectors.
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#include "support_collector.h"
#include "dirichlet/support.h"

namespace fem::bc {

SupportCollector::SupportCollector(const std::string& name)
    : model::Collection<Dirichlet::Ptr>(name) {}

constraint::Equations SupportCollector::get_equations(model::ModelData& model_data) {
    return get_equations<Support>(model_data);
}

} // namespace fem::bc
