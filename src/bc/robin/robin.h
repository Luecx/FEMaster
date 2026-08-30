/**
 * @file robin.h
 * @brief Declares Robin boundary conditions and their symbolic equation rows.
 */

#pragma once

#include "../load.h"

#include "../../core/logging.h"
#include "../../core/types_num.h"
#include "../../data/field.h"

#include <memory>
#include <vector>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

struct RobinEquationEntry {
    ID        node_id{};
    Dim       dof{};
    Precision coeff{};
};

struct RobinEquation {
    ID                              row_node_id{};
    Dim                             row_dof{};
    std::vector<RobinEquationEntry> entries{};
};

using RobinEquations = std::vector<RobinEquation>;

/**
 * @brief Boundary condition producing a RHS contribution and symbolic rows.
 */
struct Robin : Load {
    using Ptr = std::shared_ptr<Robin>;

    virtual ~Robin() = default;

    // Robin conditions must never be evaluated through the RHS-only load path.
    void apply(model::ModelData& model_data,
               model::Field&     rhs,
               Precision         time,
               bool              ignore_amplitude = false) final {
        (void) model_data;
        (void) rhs;
        (void) time;
        (void) ignore_amplitude;
        logging::error(false,
            "Robin boundary conditions require symbolic equation assembly");
    }

    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       RobinEquations&   equations,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;
};

} // namespace fem::bc
