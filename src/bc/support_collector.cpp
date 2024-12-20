/******************************************************************************
* @file support_collector.cpp
 * @brief Implements the SupportCollector class for managing collections of supports in an FEM model.
 *
 * @details The SupportCollector class provides methods for adding supports to various regions
 *          and applying them collectively to the FEM model.
 *
 * @date Created on 27.11.2024
 * @author Finn Eggers
 ******************************************************************************/

#include "support_collector.h"
#include "../data/node_data_dict.h"

namespace fem {

// Constructor
SupportCollector::SupportCollector(const std::string& p_name)
    : model::Collection<Support>(p_name, true, false) {}

// Applies all supports in the collector to the FEM model
void SupportCollector::apply(model::ModelData& model_data, NodeData& bc, constraint::Equations& equations) {
    for (Support& sup : _data) {
        sup.apply(model_data, bc, equations);
    }
}

// Adds a Support for a node region
void SupportCollector::add_supp(model::NodeRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    _data.emplace_back(Support{region, values, coordinate_system});
}

// Adds a Support for an element region
void SupportCollector::add_supp(model::ElementRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    _data.emplace_back(Support{region, values, coordinate_system});
}

// Adds a Support for a surface region
void SupportCollector::add_supp(model::SurfaceRegion::Ptr region, Vec6 values, cos::CoordinateSystem::Ptr coordinate_system) {
    _data.emplace_back(Support{region, values, coordinate_system});
}

} // namespace fem
