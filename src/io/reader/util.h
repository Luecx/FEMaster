/**
 * @file util.h
 * @brief Declares shared reader utilities for nodal boundary-condition regions.
 *
 * The reader resolves command-wide and nodal coordinate-system assignments
 * before creating loads or supports. These utilities group compiled nodes by
 * their effective basis; reference resolution remains in model::Model and basis
 * evaluation during assembly remains in the boundary-condition implementation.
 *
 * @see Parser
 * @see group_by_orientation
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#pragma once

#include "../../cos/coordinate_system.h"
#include "../../data/region.h"

#include <utility>
#include <vector>

namespace fem::io::reader {

class Parser;

// Group compiled nodes by nodal TRANSFORM or the supplied default basis.
// Uniform regions are preserved; a null orientation denotes the global basis.
std::vector<std::pair<model::NodeRegion::Ptr, cos::CoordinateSystem::Ptr>>
group_by_orientation(const Parser&              parser,
                     model::NodeRegion::Ptr      region,
                     cos::CoordinateSystem::Ptr default_orientation);

} // namespace fem::io::reader
