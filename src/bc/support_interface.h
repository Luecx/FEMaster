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
 * @see dirichlet/dirichlet.h
 * @see Support
 * @see Temperature
 * @see SupportCollector
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "dirichlet/dirichlet.h"

#include "../constraints/types/equation.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Solver-equation interface shared by Dirichlet support conditions.
 *
 * Mechanical and thermal definitions may coexist in one collector even though
 * local degree of freedom zero denotes x-translation in the structural system
 * and temperature in the thermal system. The requesting load case therefore
 * selects entries by their concrete type.
 */
struct SupportInterface : Dirichlet {
    using Ptr = std::shared_ptr<SupportInterface>;

    ~SupportInterface() override = default;

    virtual void apply(model::ModelData& model_data,
                       constraint::Equations& equations) = 0;

    std::string str() const override = 0;
};

} // namespace fem::bc
