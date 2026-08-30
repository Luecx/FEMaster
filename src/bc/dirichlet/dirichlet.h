/**
 * @file dirichlet.h
 * @brief Declares the common interface for essential boundary conditions.
 *
 * Dirichlet conditions prescribe primary solution variables and therefore enter
 * the algebraic system as constraint equations rather than right-hand-side
 * contributions. Mechanical supports and prescribed temperatures share this
 * contract while retaining solver-specific concrete representations.
 *
 * @see Support
 * @see Temperature
 * @see SupportCollector
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "../bc.h"

#include "../../constraints/types/equation.h"
#include "../../core/printable.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common abstraction for prescribed primary solution variables.
 *
 * A Dirichlet boundary condition translates its physical prescription into one
 * or more equations over compiled model entities. The requesting load case owns
 * the interpretation of the local degree of freedom: structural supports use
 * generalized displacement components, whereas thermal temperature conditions
 * use scalar thermal component zero.
 *
 * Concrete conditions own their target region and prescribed values. They append
 * complete equations in deterministic order and provide a compact diagnostic
 * representation through `Printable`.
 */
struct Dirichlet : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Dirichlet>;

    virtual ~Dirichlet() = default;

    // Translate the prescribed primary variable into algebraic equations for
    // the compiled model. Concrete conditions append to the supplied container.
    virtual void apply(model::ModelData& model_data, constraint::Equations& equations) = 0;

    // Return the concrete condition in a compact diagnostic representation.
    std::string str() const override = 0;
};

} // namespace fem::bc
