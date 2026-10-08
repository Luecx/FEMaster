/**
 * @file util.cpp
 * @brief Implements grouping of compiled nodes by their effective load basis.
 *
 * Nodal TRANSFORM assignments stored by Parser override command-wide defaults.
 * Grouping uses coordinate-system pointer identity, preserves the original
 * region when its basis is uniform and otherwise creates private node regions.
 * The model owns coordinate-system definitions; this utility neither transforms
 * components nor modifies parser state or the input region.
 *
 * @see group_by_orientation
 * @see Parser
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#include "util.h"

#include "parser.h"
#include "../../core/logging.h"
#include "../../model/model.h"

#include <memory>

namespace fem::io::reader {

/**
 * @brief Groups a compiled node region by its effective coordinate system.
 *
 * Every node uses its nodal *TRANSFORM assignment when one exists. Otherwise,
 * the supplied default orientation is used, which may be null to represent the
 * global coordinate system. Nodes with the same effective coordinate system are
 * collected into one region.
 *
 * The input region itself is returned unchanged when every node uses the same
 * effective orientation. This preserves existing named regions in the common
 * case where no splitting is required. An empty region produces no groups;
 * a null region or an unknown assigned coordinate system is rejected.
 *
 * Groups follow the first occurrence of their orientation in the input region.
 * Coordinate systems are compared by identity, so distinct definitions remain
 * separate even if their bases are geometrically equivalent. Parser state and
 * the input region are not modified.
 *
 * @param parser              Parser providing the model and nodal assignments.
 * @param region              Compiled node region to group.
 * @param default_orientation Coordinate system used for nodes without *TRANSFORM.
 * @return Node regions paired with their effective coordinate systems.
 */
std::vector<std::pair<model::NodeRegion::Ptr, cos::CoordinateSystem::Ptr>>
group_by_orientation(const Parser&              parser,
                     model::NodeRegion::Ptr      region,
                     cos::CoordinateSystem::Ptr default_orientation) {
    // Validate the region before accessing its nodes or allocating local state.
    logging::error(region != nullptr,
        "group_by_orientation: node region is not initialized");

    const auto& model = parser.model();

    // Resolve the effective orientation of every node. Nodal *TRANSFORM
    // assignments take precedence over the orientation supplied by the command.
    std::vector<std::pair<ID, cos::CoordinateSystem::Ptr>> node_orientations;
    node_orientations.reserve(region->size());

    for (const ID node_id : *region) {
        cos::CoordinateSystem::Ptr orientation = default_orientation;

        const auto transform = parser.node_transforms.find(node_id);
        if (transform != parser.node_transforms.end()) {
            logging::error(model._data->coordinate_systems.has(transform->second),
                "Node ", node_id, " references unknown coordinate system ", transform->second);

            orientation = model._data->coordinate_systems.get(transform->second);
        }

        node_orientations.emplace_back(node_id, std::move(orientation));
    }

    // Preserve the original region when all nodes use the same orientation.
    // This avoids constructing private regions for the common unsplit case.
    if (!node_orientations.empty()) {
        const auto& orientation = node_orientations.front().second;

        bool same_orientation = true;
        for (const auto& node_orientation : node_orientations) {
            if (node_orientation.second != orientation) {
                same_orientation = false;
                break;
            }
        }

        if (same_orientation) {
            return {{std::move(region), orientation}};
        }
    }

    // Group nodes by coordinate-system identity. Different TRANSFORM definitions
    // remain separate even if their bases happen to be geometrically equivalent.
    std::vector<std::pair<model::NodeRegion::Ptr, cos::CoordinateSystem::Ptr>> groups;

    for (const auto& [node_id, orientation] : node_orientations) {
        auto group = groups.end();

        for (auto it = groups.begin(); it != groups.end(); ++it) {
            if (it->second == orientation) {
                group = it;
                break;
            }
        }

        if (group == groups.end()) {
            auto grouped_region = std::make_shared<model::NodeRegion>("INTERNAL");
            grouped_region->add(node_id);

            groups.emplace_back(std::move(grouped_region), orientation);
        } else {
            group->first->add(node_id);
        }
    }

    return groups;
}

} // namespace fem::io::reader
