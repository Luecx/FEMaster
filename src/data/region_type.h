/**
 * @file region_type.h
 * @brief Defines compile-time categories for model identifier regions.
 *
 * The model-data subsystem tags Region specializations with node, element,
 * surface or line identity. The category determines which model entities an
 * identifier collection denotes; it neither stores entities nor validates IDs.
 * Collection and Region implement storage and diagnostics separately.
 *
 * @see RegionTypes
 * @see Region
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#pragma once

namespace fem {
namespace model {

// Underlying representation for region-kind tags and numeric diagnostics.
using RegionType = int;

/**
 * @brief Distinguishes the model entity kinds represented by typed regions.
 *
 * These tags are template arguments of Region and are independent of FieldDomain,
 * which describes numerical field-row layouts rather than identifier collections.
 */
enum RegionTypes : RegionType {
    // Node identifiers.
    NODE,
    // Element identifiers.
    ELEMENT,
    // Surface identifiers.
    SURFACE,
    // Line identifiers for one-dimensional geometry.
    LINE
};
} // namespace model
} // namespace fem
