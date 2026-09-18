/**
 * @file thermal.h
 * @brief Defines the physical-domain marker shared by thermal boundary conditions.
 *
 * FEMaster classifies boundary conditions algebraically through Dirichlet,
 * Neumann and Mixed. Thermal conditions cut across those categories:
 * prescribed temperature is Dirichlet, prescribed heat flux is Neumann and
 * convection is Mixed.
 *
 * `ThermalCondition` therefore does not introduce a fourth algebraic boundary-
 * condition category and deliberately does not derive from `BoundaryCondition`.
 * It is an orthogonal polymorphic marker used only to guarantee that a
 * `ThermalCollector` contains thermal definitions.
 *
 * @see Dirichlet
 * @see Neumann
 * @see Mixed
 * @see ThermalCollector
 *
 * @author Finn Eggers
 * @date 18.09.2026
 */

#pragma once

#include <memory>

namespace fem::bc {

/**
 * @brief Polymorphic marker identifying boundary conditions of the thermal field.
 *
 * A concrete thermal condition also derives from exactly one algebraic boundary-
 * condition category. The two inheritance axes describe independent concepts:
 * the category determines how the condition enters the discrete system, while
 * this marker determines that the primary physical field is temperature.
 *
 * The marker intentionally exposes no assembly operation. `ThermalCollector`
 * dispatches each stored object through its Dirichlet, Load or Mixed interface.
 */
struct ThermalCondition {
    // Types
    using Ptr = std::shared_ptr<ThermalCondition>;

    // Polymorphic destruction required for cross-casting to algebraic BC types
    virtual ~ThermalCondition() = default;
};

} // namespace fem::bc
