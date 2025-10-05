/**
 * @file region_type.h
 * @brief Declares the enumeration that categorises region collections.
 *
 * Region types are used to distinguish node, element, and surface regions when
 * parsing and assembling model data.
 *
 * @see src/data/region_type.cpp
 * @see src/data/region.h
 */

#pragma once

namespace fem {
namespace model {

using RegionType = int; ///< Underlying type used by the `RegionTypes` enumeration.

/// Enumerates supported region kinds recognised by the model.
enum RegionTypes : RegionType {
    NODE,    ///< Collection of node identifiers.
    ELEMENT, ///< Collection of element identifiers.
    SURFACE  ///< Collection of surface identifiers.
};

} // namespace model
} // namespace fem

