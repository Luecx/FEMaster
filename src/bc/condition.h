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
 *     f   : nodal right-hand-side field,
 *     C,d : algebraic constraint equations,
 *     K_c : condition-dependent matrix contribution.
 *
 * A concrete condition modifies only the outputs required by its physical
 * definition and leaves the remaining outputs unchanged. This keeps the common
 * interface general without introducing separate Dirichlet, Neumann, Mixed or
 * load-capability hierarchies.
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

    // Polymorphic lifetime shared between active storage and named collectors
    ~Condition() override = default;

    // Assemble into rhs, constraints or matrix using the caller's multiplier.
    virtual void apply(
        model::ModelData&      model_data,
        model::Field&          rhs,
        constraint::Equations& equations,
        const SystemDofIds&    system_dof_ids,
        TripletList&           matrix,
        Precision              scale = Precision(1)
    ) = 0;

    // Diagnostics
    std::string str() const override = 0;
};

} // namespace fem::bc
