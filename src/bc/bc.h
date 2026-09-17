/**
 * @file bc.h
 * @brief Defines the common polymorphic root of FEMaster boundary conditions.
 *
 * The boundary-condition subsystem separates prescribed kinematic or thermal
 * values, externally applied loads and conditions that contribute to both the
 * right-hand side and the system operator. The concrete mathematical behavior
 * is provided by the `Dirichlet`, `Neumann` and `Mixed` category bases in their
 * respective subdirectories.
 *
 * `BoundaryCondition` itself deliberately carries no assembly interface. It is
 * the semantic root used to identify all boundary-condition definitions while
 * the category bases expose the operations appropriate to their algebraic role.
 *
 * @see Dirichlet
 * @see Neumann
 * @see Mixed
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include <memory>

namespace fem::bc {

/**
 * @brief Semantic root of all FEMaster boundary-condition definitions.
 *
 * A boundary condition modifies the unconstrained finite-element problem in one
 * of three ways:
 *
 * - Dirichlet conditions prescribe primary variables through algebraic
 *   constraint equations,
 * - Neumann conditions contribute only to the assembled right-hand side,
 * - mixed conditions contribute to both the right-hand side and the system
 *   operator.
 *
 * This base class intentionally defines only polymorphic ownership and
 * destruction. Derived category bases provide the actual assembly contracts so
 * code cannot accidentally treat every boundary condition as the same algebraic
 * operation.
 */
struct BoundaryCondition {
    // Types
    using Ptr = std::shared_ptr<BoundaryCondition>;

    // Polymorphic destruction
    virtual ~BoundaryCondition() = default;
};

} // namespace fem::bc
