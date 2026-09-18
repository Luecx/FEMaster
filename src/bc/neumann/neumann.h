/**
 * @file neumann.h
 * @brief Defines the common base for right-hand-side-only boundary conditions.
 *
 * FEMaster uses the Neumann category for conditions whose discrete contribution
 * enters only the generalized external load vector. Classical surface tractions
 * are the canonical Neumann boundary term, but the same algebraic category also
 * contains volumetric, inertial and equivalent thermal loads because they share
 * the same solver-facing assembly contract.
 *
 * In a weak finite-element equilibrium statement
 *
 *     delta u^T f_int(u) = delta u^T f_ext,
 *
 * a Neumann condition contributes only to `f_ext`; it does not prescribe a
 * primary variable and does not add an independent system-operator term.
 *
 * @see BoundaryCondition
 * @see Load
 * @see Mixed
 * @see HeatFlux
 *
 * @author Finn Eggers
 * @date 17.09.2026
 */

#pragma once

#include "../bc.h"
#include "../load.h"

#include <memory>

namespace fem::bc {

/**
 * @brief Base class for conditions contributing only to the external RHS.
 *
 * A concrete Neumann condition converts its physical definition to equivalent
 * generalized nodal loads through `Load::apply()`. Depending on the condition,
 * this conversion may be a direct nodal assignment or a consistent finite-
 * element integral such as
 *
 *     f_e = integral_Gamma N^T t dGamma
 *
 * for a surface traction, or
 *
 *     f_e = integral_Omega N^T b dOmega
 *
 * for a distributed body-force density.
 *
 * The category intentionally describes the algebraic solver contribution rather
 * than only the strict PDE boundary classification. Consequently, body,
 * inertia and equivalent thermal loads remain Neumann conditions in the
 * FEMaster source architecture because they modify only the assembled RHS.
 */
struct Neumann : BoundaryCondition, Load {
    // Types
    using Ptr = std::shared_ptr<Neumann>;

    // Polymorphic destruction
    ~Neumann() override = default;
};

} // namespace fem::bc
