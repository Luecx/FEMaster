/******************************************************************************
 * @file config.h
 * @brief Declares global configuration parameters for the FEM runtime.
 *
 * Provides a lightweight structure that stores process-wide settings such as
 * the number of worker threads and exposes a singleton-style instance.
 *
 * @see src/core/config.cpp
 * @see src/core/types_num.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#pragma once

#include "types_num.h"

namespace fem {

/******************************************************************************
 * @struct Config
 * @brief Holds shared runtime configuration values.
 ******************************************************************************/
struct Config {
    ID max_threads = 1; ///< Maximum number of threads used by the solver.
};

/******************************************************************************
 * @brief Global configuration object accessible throughout the codebase.
 ******************************************************************************/
extern Config global_config;

} // namespace fem
