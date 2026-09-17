/**
 * @file mixed.h
 * @brief Defines the common base for boundary conditions with RHS and operator terms.
 *
 * Mixed conditions couple the primary variable to a boundary flux or load and
 * therefore contribute to both sides of the discrete finite-element system. A
 * classical example is a Robin law such as convection,
 *
 *     q_n = h (T - T_inf),
 *
 * whose weak form produces a boundary matrix proportional to
 *
 *     integral_Gamma h N^T N dGamma
 *
 * and a prescribed right-hand-side term proportional to
 *
 *     integral_Gamma h T_inf N^T dGamma.
 *
 * FEMaster keeps these two contributions explicit: `Load::apply()` assembles the
 * RHS part and `Mixed::apply_matrix()` assembles the sparse operator part.
 * Constraint equations are not involved.
 *
 * @see BoundaryCondition
 * @see Load
 * @see Neumann
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
 * @brief Base class for conditions contributing to both RHS and system operator.
 *
 * A mixed condition represents a constitutive relation on a boundary rather
 * than a prescribed primary variable or a prescribed flux alone. After spatial
 * discretization, the condition generally contributes terms of the form
 *
 *     K_b u = f_b,
 *
 * where `K_b` is assembled from the part depending on the unknown field and
 * `f_b` from the prescribed part. The surrounding load-case assembly determines
 * the final sign with which these terms enter the global residual.
 *
 * Derived classes implement the normal `Load::apply()` contract for `f_b` and
 * additionally provide sparse triplets for `K_b`. The `system_dof_ids` mapping
 * converts model-level nodal DOFs to the active system numbering used by the
 * global matrix.
 */
struct Mixed : BoundaryCondition, Load {
    // Types
    using Ptr = std::shared_ptr<Mixed>;

    // Polymorphic destruction
    ~Mixed() override = default;

    // Assemble the operator part of the mixed condition into sparse system
    // triplets. RHS assembly remains provided by the inherited Load interface.
    virtual void apply_matrix(model::ModelData&     model_data,
                              const SystemDofIds&   system_dof_ids,
                              TripletList&          matrix,
                              Precision             time,
                              bool                  ignore_amplitude = false) = 0;
};

} // namespace fem::bc
