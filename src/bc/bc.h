/**
 * @file bc.h
 * @brief Declares the common root type for boundary-condition objects.
 *
 * The boundary-condition module contains loads, supports and their associated
 * collectors. `BoundaryCondition` provides the minimal polymorphic ownership
 * type shared by these objects without imposing assembly semantics on derived
 * classes. More specialized interfaces, such as `Load`, add the operations
 * required by their respective solver pipelines.
 *
 * @see BoundaryCondition
 * @see load.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <memory>

namespace fem {
namespace bc {

/**
 * @brief Serves as the common ownership root of boundary-condition types.
 *
 * The type deliberately contains no state and currently exposes only a shared
 * pointer alias. It establishes a stable extension point for functionality that
 * is common to every boundary condition while allowing specialized classes to
 * define their own application interfaces.
 */
struct BoundaryCondition {
    // Shared ownership type used when boundary conditions are stored
    // polymorphically by model and load-case data structures.
    using Ptr = std::shared_ptr<BoundaryCondition>;
};
} // namespace bc
} // namespace fem
