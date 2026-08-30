/**
 * @file load.h
 * @brief Declares the common ownership base for load-side boundary conditions.
 */

#pragma once

#include "amplitude.h"
#include "bc.h"

#include "../core/printable.h"
#include "../cos/coordinate_system.h"

#include <memory>
#include <string>

namespace fem::bc {

/**
 * @brief Common state shared by Neumann and Robin boundary conditions.
 *
 * `Load` intentionally has no assembly method. Neumann conditions contribute
 * only to the right-hand side, whereas Robin conditions additionally emit
 * symbolic linear equation rows. Their distinct interfaces live in the
 * corresponding derived base classes.
 */
struct Load : BoundaryCondition, Printable {
    using Ptr = std::shared_ptr<Load>;

    cos::CoordinateSystem::Ptr orientation_ = nullptr;
    Amplitude::Ptr amplitude_ = nullptr;

    virtual ~Load() = default;
    std::string str() const override = 0;
};

} // namespace fem::bc
