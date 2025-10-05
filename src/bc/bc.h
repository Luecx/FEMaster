/******************************************************************************
 * @file bc.h
 * @brief Declares the common base type for all boundary conditions.
 *
 * The `BoundaryCondition` struct provides a shared interface for all boundary
 * condition types in the boundary-condition (BC) module. It is primarily a
 * symbolic base used for polymorphic ownership and future extensibility.
 *
 * @see src/bc/load.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include <memory>

namespace fem {
namespace bc {

/******************************************************************************
 * @struct BoundaryCondition
 * @brief Serves as the abstract root for BC polymorphism.
 *
 * This base struct currently stores only the shared-pointer alias that is used
 * by collectors. Concrete boundary condition implementations derive from this
 * type.
 ******************************************************************************/
struct BoundaryCondition {
    using Ptr = std::shared_ptr<BoundaryCondition>; ///< Shared pointer alias for BC ownership.
};

} // namespace bc
} // namespace fem
