/**
 * @file support_collector.cpp
 * @brief Implements aggregation of support boundary conditions.
 *
 * The collector offers convenience routines for defining supports over various
 * regions and producing constraint equations in a single call.
 *
 * @see src/bc/support_collector.h
 * @see src/bc/support.cpp
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "support_collector.h"

#include "../data/node_data_dict.h"

#include <utility>

namespace fem {
namespace bc {

/**
 * @copydoc SupportCollector::SupportCollector
 */
SupportCollector::SupportCollector(const std::string& name)
    : model::Collection<Support>(name, true, false) {}

/**
 * @copydoc SupportCollector::get_equations
 */
constraint::Equations SupportCollector::get_equations(model::ModelData& model_data) {
    constraint::Equations equations{};
    for (Support& support : this->_data) {
        support.apply(model_data, equations);
    }
    return equations;
}

/**
 * @copydoc SupportCollector::add_supp(model::NodeRegion::Ptr,Vec6,cos::CoordinateSystem::Ptr)
 */
void SupportCollector::add_supp(model::NodeRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    this->_data.emplace_back(Support{std::move(region), values, std::move(coordinate_system)});
}

/**
 * @copydoc SupportCollector::add_supp(model::ElementRegion::Ptr,Vec6,cos::CoordinateSystem::Ptr)
 */
void SupportCollector::add_supp(model::ElementRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    this->_data.emplace_back(Support{std::move(region), values, std::move(coordinate_system)});
}

/**
 * @copydoc SupportCollector::add_supp(model::SurfaceRegion::Ptr,Vec6,cos::CoordinateSystem::Ptr)
 */
void SupportCollector::add_supp(model::SurfaceRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    this->_data.emplace_back(Support{std::move(region), values, std::move(coordinate_system)});
}

} // namespace bc
} // namespace fem
