/**
 * @file neumann.h
 * @brief Declares the common base for natural boundary conditions.
 */

#pragma once

#include "../bc.h"
#include "../../core/printable.h"

#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Common ownership and diagnostics interface for Neumann conditions.
 *
 * Neumann conditions contribute generalized fluxes to a residual and may also
 * contribute boundary tangent terms, as for thermal convection.
 */
struct Neumann : BoundaryCondition, fem::Printable {
    using Ptr = std::shared_ptr<Neumann>;

    ~Neumann() override = default;
    std::string str() const override = 0;
};

} // namespace fem::bc
