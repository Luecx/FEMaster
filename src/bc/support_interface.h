/**
 * @file support_interface.h
 * @brief Declares the common interface for prescribed solution variables.
 *
 * Support conditions impose essential boundary values on one physical solver
 * system. Mechanical supports prescribe generalized displacements, whereas
 * thermal supports prescribe scalar nodal temperatures. Both kinds share named
 * support collectors, while each load case selects the concrete support type
 * compatible with its own degree-of-freedom mapping.
 *
 * The interface deliberately owns no region or prescribed value. Concrete
 * support types retain those definitions and translate them into algebraic
 * equations compatible with their solver-specific degree-of-freedom mapping.
 *
 * @see Support
 * @see Temperature
 * @see SupportCollector
 *
 * @author Finn Eggers
 * @date 29.08.2026
 */

#pragma once

#include "bc.h"

#include "../constraints/types/equation.h"
#include "../core/printable.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common polymorphic contract for essential boundary conditions.
 *
 * A support interface converts one prescribed primary variable into constraint
 * equations over compiled model nodes. Mechanical and thermal definitions may
 * coexist in one collector even though local degree of freedom zero denotes
 * x-translation in the structural system and temperature in the thermal system.
 * The requesting load case therefore selects entries by their concrete type.
 *
 * Derived classes own their target region, prescribed values and any coordinate
 * conventions. They must append complete equations in deterministic order and
 * provide a concise printable representation for diagnostics.
 */
struct SupportInterface : BoundaryCondition, Printable {
    // Shared ownership type used by heterogeneous support collectors
    using Ptr = std::shared_ptr<SupportInterface>;

    // Polymorphic destruction through collector-owned base pointers
    virtual ~SupportInterface() = default;

    // Translate the concrete prescribed values into solver equations over the
    // compiled model topology.
    virtual void apply(model::ModelData& model_data, constraint::Equations& equations) = 0;

    // Return the concrete support definition for diagnostics
    std::string str() const override = 0;
};

} // namespace fem::bc
