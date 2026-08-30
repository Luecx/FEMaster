/**
 * @file dirichlet.h
 * @brief Declares the common base for essential boundary conditions.
 */

#pragma once

#include "../bc.h"
#include "../../core/printable.h"

#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Common ownership and diagnostics interface for Dirichlet conditions.
 *
 * Dirichlet conditions prescribe primary solution variables. Concrete solver
 * domains provide the equation construction required for their unknowns.
 */
struct Dirichlet : BoundaryCondition, fem::Printable {
    using Ptr = std::shared_ptr<Dirichlet>;

    ~Dirichlet() override = default;
    std::string str() const override = 0;
};

} // namespace fem::bc
