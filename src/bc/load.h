/**
 * @file load.h
 * @brief Declares the common ownership base for load-side boundary conditions.
 */

#pragma once

#include "amplitude.h"
#include "bc.h"

#include "../core/printable.h"
#include "../core/types_num.h"
#include "../cos/coordinate_system.h"
#include "../data/field.h"

#include <memory>
#include <string>

namespace fem::model {
struct ModelData;
}

namespace fem::bc {

/**
 * @brief Common state and RHS dispatch shared by load-side conditions.
 *
 * Neumann conditions implement the RHS dispatch directly. Robin conditions
 * deliberately reject that generic path and must be evaluated through their
 * equation-producing overload.
 */
struct Load : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Load>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr amplitude_ = nullptr;

    virtual ~Load() = default;

    virtual void apply(model::ModelData& model_data,
                       model::Field&     rhs,
                       Precision         time,
                       bool              ignore_amplitude = false) = 0;

    std::string str() const override = 0;
};

} // namespace fem::bc
