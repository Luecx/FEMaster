/**
 * @file condition.h
 * @brief Defines the common interface for FEMaster model conditions.
 *
 * A condition describes externally prescribed behavior in the currently
 * active model state. Structural supports, structural loads and thermal boundary
 * conditions all share this semantic layer even though they contribute to
 * different parts of the discrete finite-element system.
 *
 * Every condition receives the complete set of solver-facing assembly outputs:
 *
 *     f     : nodal external load field,
 *     K_L   : optional load tangent df/du,
 *     C,d   : algebraic constraint equations,
 *     K_c   : condition-dependent system-matrix contribution.
 *
 * A concrete condition modifies only the outputs required by its physical
 * definition and leaves the remaining outputs unchanged. This keeps the common
 * interface general without introducing separate Dirichlet, Neumann, Mixed or
 * load-capability hierarchies.
 *
 * The system matrix K_c and the load tangent K_L are deliberately separate.
 * K_c is part of the system operator in every analysis, for example the
 * convection film operator. K_L is the analytic derivative of a follower load
 * with respect to the nodal displacements in the current configuration and
 * only enters the Newton linearization of a nonlinear analysis:
 *
 *     K_eff = K_t + K_c - K_L.
 *
 * A condition always writes the true derivative df/du into K_L; the calling
 * analysis subtracts it. Dead loads are evaluated in the reference
 * configuration and never contribute to K_L. A follower load requires a load
 * tangent and reports an error when the caller does not request one, because
 * its residual would otherwise be linearized inconsistently.
 *
 * Condition history remains independent of assembly. The load case owns independent initial and target condition states, while
 * the reader resolves original identifiers and semantic replacement rules.
 * Named collectors retain reusable definitions independently of the active set.
 *
 * @see ConditionState
 *
 * @author Finn Eggers
 * @date 07.10.2026
 */

#pragma once

#include "amplitude.h"

#include "../constraints/types/equation.h"
#include "../core/printable.h"
#include "../core/types_eig.h"
#include "../core/types_num.h"
#include "../model/system_dofs.h"

#include <cmath>
#include <memory>
#include <string>

namespace fem::model {
struct Field;
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Physical boundary-condition or load definition.
 *
 * LoadCase owns the initial and target states; the parser handles replacement.
 * Physical regions, magnitudes and optional amplitude references stay on the
 * definition. The analysis supplies the numerical scale for additive loads;
 * absolute Dirichlet conditions ignore it. Time and interpolation do not belong
 * in a physical condition implementation.
 */
struct Condition : Printable {
    // Shared ownership retains physical definitions in model history or reusable
    // named groups without coupling their lifetime to an analysis object.
    using Ptr = std::shared_ptr<Condition>;

    // Optional reusable amplitude sampled by the calling assembly.
    Amplitude::Ptr amplitude_ = nullptr;

    // Configuration in which the condition is evaluated. A follower load reads
    // the current positions and contributes an analytic load tangent; a dead
    // load reads the reference positions and is independent of the displacement.
    bool follower_ = false;

    // Polymorphic lifetime shared between active storage and named collectors
    ~Condition() override = default;

    // Assemble into the external load, its optional tangent, the constraints or
    // the system matrix using the caller's multiplier. Load tangent triplets use
    // the global equation indices of system_dof_ids and contain df/du with the
    // same multiplier as the load. A null load_tangent requests no linearization.
    virtual void apply(
        model::ModelData&        model_data,
        model::Field&            load,
        TripletList*             load_tangent,
        constraint::Equations&   constraints,
        const model::SystemDofs& system_dof_ids,
        TripletList&             system_matrix,
        Precision                scale
    ) = 0;

    // Diagnostics
    std::string str() const override = 0;

protected:
    // Adds the 3x3 load-tangent block df(row_node, row_dof + i) / dx(col_node, j)
    // to the triplets. Nodal positions only change through the translational
    // DOFs, so the columns always address components 0..2. Components without
    // an active equation do not exist in the system and are skipped.
    static void add_load_tangent_block(TripletList&             load_tangent,
                                       const model::SystemDofs& system_dof_ids,
                                       Index                    row_node,
                                       Index                    row_dof,
                                       Index                    col_node,
                                       const Mat3&              block) {
        for (Index i = 0; i < 3; ++i) {
            const Index row = system_dof_ids(row_node, row_dof + i);
            if (row < 0) continue;

            for (Index j = 0; j < 3; ++j) {
                const Index col = system_dof_ids(col_node, j);
                if (col < 0 || block(i, j) == Precision(0)) continue;

                load_tangent.emplace_back(row, col, block(i, j));
            }
        }
    }
};

} // namespace fem::bc
