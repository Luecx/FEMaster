/**
 * @file load.h
 * @brief Declares common state for load-side boundary conditions.
 */

#pragma once

#include "amplitude.h"
#include "bc.h"

#include "../core/printable.h"
#include "../core/types_cls.h"
#include "../cos/coordinate_system.h"

#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Common ownership and modifier state for load-side conditions.
 *
 * The base deliberately defines no assembly operation. Neumann and Robin
 * conditions have different algebraic contracts and expose their own `apply()`
 * interfaces in the corresponding derived base classes.
 */
struct Load : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Load>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr amplitude_ = nullptr;

    virtual ~Load() = default;
    std::string str() const override = 0;
};

} // namespace fem::bc
