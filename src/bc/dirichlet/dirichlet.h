/**
 * @file dirichlet.h
 * @brief Defines the common base for prescribed primary-variable boundary conditions.
 *
 * Dirichlet conditions prescribe finite-element primary variables instead of
 * adding physical load terms. FEMaster represents these prescriptions as
 * algebraic equations collected independently from the unconstrained system and
 * consumed by the constraint-transformation machinery.
 *
 * For a discrete unknown vector `u`, a set of Dirichlet conditions is written
 * as
 *
 *     C u = d,
 *
 * where every generated equation contributes one row of the constraint matrix
 * `C` and one entry of the prescribed-value vector `d`. A direct global support
 * produces a unit row, whereas an oriented support may couple several global
 * components through the local basis directions.
 *
 * @see BoundaryCondition
 * @see constraint::Equation
 * @see Support
 * @see Temperature
 *
 * @author Finn Eggers
 * @date 17.09.2026
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
 * @brief Base class for essential conditions represented by constraint equations.
 *
 * A Dirichlet condition does not directly modify the finite-element stiffness
 * matrix or external load vector. Instead it appends equations of the form
 *
 *     sum_j C_ij u_j = d_i
 *
 * to a `constraint::Equations` collection. The global constraint handling then
 * enforces these rows together with other kinematic relations when the active
 * system is formed.
 *
 * Derived classes are responsible for resolving their semantic target to global
 * degrees of freedom and for constructing the correct coefficients in `C`.
 * Structural supports may therefore generate rows in the six generalized nodal
 * DOFs, while prescribed temperature conditions use the scalar thermal DOF.
 */
struct Dirichlet : BoundaryCondition, Printable {
    // Types
    using Ptr = std::shared_ptr<Dirichlet>;

    // Polymorphic destruction
    ~Dirichlet() override = default;

    // Resolve the condition target and append its algebraic rows C_i u = d_i
    // to the supplied equation collection.
    virtual void apply(model::ModelData&       model_data,
                       constraint::Equations& equations) = 0;

    // Diagnostics
    std::string str() const override = 0;
};

} // namespace fem::bc
