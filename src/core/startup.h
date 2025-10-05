/**
 * @file startup.h
 * @brief Declares the startup helper that prints build information at launch.
 *
 * Instantiating the `Startup` struct emits version and capability data to the
 * log output.
 *
 * @see src/core/startup.cpp
 * @see src/core/version.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

namespace fem {
namespace startup {

/**
 * @struct Startup
 * @brief Emits build information during static initialisation.
 */
struct Startup {
    Startup();
};

extern Startup instance; ///< Global startup object that triggers automatic reporting.

} // namespace startup
} // namespace fem
