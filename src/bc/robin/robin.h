/**
 * @file robin.h
 * @brief Declares mixed boundary conditions contributing to both RHS and LHS.
 *
 * Robin conditions contain a prescribed source part and a part proportional to
 * the primary solution variable on the boundary. Their weak finite-element form
 * therefore contributes simultaneously to the external right-hand side and to
 * the global operator. The current implementation assembles operator terms as
 * sparse triplets using an analysis-specific scalar DOF map.
 *
 * The hierarchy derives from `Load` so Robin and Neumann conditions can share
 * collector ownership, amplitudes and common diagnostics. The RHS-only `Load`
 * dispatch is deliberately rejected to prevent a caller from silently dropping
 * the operator contribution of a mixed boundary condition.
 *
 * @see Load
 * @see Robin
 * @see ../neumann/neumann.h
 * @see convection.h
 *
 * @author Finn Eggers
 * @date 30.08.2026
 */

#pragma once

#include "../load.h"

#include "../../core/logging.h"
#include "../../core/types_eig.h"
#include "../../core/types_num.h"
#include "../../data/field.h"

#include <memory>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common abstraction for boundary conditions contributing to RHS and LHS.
 *
 * A Robin condition has the generic linear boundary form
 *
 * \f[
 *     q = g - \beta u,
 * \f]
 *
 * where `u` is the primary boundary unknown. In the weak formulation, `g`
 * generates a consistent nodal source vector while the `beta u` term generates
 * a boundary matrix. Concrete implementations own the physical field and local
 * surface integration; the caller supplies the active-system DOF mapping needed
 * to convert local nodal matrix entries into sparse global triplets.
 *
 * The class intentionally finalizes the inherited RHS-only `apply()` overload
 * with an error. This makes it impossible to treat a Robin condition as an
 * ordinary Neumann load and accidentally assemble only half of its weak form.
 */
struct Robin : Load {
    // Shared ownership type used by load collectors and model-side mixed-boundary
    // dispatch.
    using Ptr = std::shared_ptr<Robin>;

    // Enable destruction through the Robin interface while retaining ownership
    // of concrete mixed boundary conditions in shared-pointer collectors.
    virtual ~Robin() = default;

    // Reject the generic RHS-only load path because a Robin condition is not
    // mathematically complete without its operator contribution.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) final {
        // Mark all common-interface arguments as intentionally unused before the
        // explicit failure. No partial RHS contribution is assembled here.
        (void) model_data;
        (void) rhs;
        (void) time;
        (void) ignore_amplitude;

        logging::error(false,
            "Robin boundary conditions require LHS assembly context");
    }

    // Assemble both parts of the mixed boundary condition. `rhs` receives the
    // prescribed source contribution; `lhs` receives global operator triplets
    // mapped through the supplied analysis-specific system DOF identifiers.
    virtual void apply(model::ModelData&   model_data,
                       model::Field&       rhs,
                       const SystemDofIds& system_dof_ids,
                       TripletList&         lhs,
                       Precision           time,
                       bool                ignore_amplitude = false) = 0;
};

} // namespace fem::bc
