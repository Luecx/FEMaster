/**
 * @file robin.h
 * @brief Declares boundary conditions contributing to both RHS and LHS.
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
 * @brief Boundary condition contributing to both the load vector and operator.
 *
 * Robin conditions assemble their prescribed source contribution into the nodal
 * right-hand side and append their state-dependent operator contribution directly
 * to a sparse triplet list. The concrete condition owns the local boundary
 * integration, while the supplied system DOF map provides the local-to-global
 * algebraic mapping.
 */
struct Robin : Load {
    using Ptr = std::shared_ptr<Robin>;

    virtual ~Robin() = default;

    // Prevent accidental use through structural RHS-only load paths.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) final {
        (void) model_data;
        (void) rhs;
        (void) time;
        (void) ignore_amplitude;
        logging::error(false,
            "Robin boundary conditions require LHS assembly context");
    }

    virtual void apply(model::ModelData&  model_data,
                       model::Field&      rhs,
                       const SystemDofIds& system_dof_ids,
                       TripletList&        lhs,
                       Precision           time,
                       bool                ignore_amplitude = false) = 0;
};

} // namespace fem::bc
