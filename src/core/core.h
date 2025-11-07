/**
 * @file core.h
 * @brief Convenience umbrella header for frequently used core utilities.
 *
 * Pulls in assertion helpers, logging, timers, and Eigen type definitions to
 * simplify translation units that rely on the core subsystem.
 *
 * @see src/core/assert.h
 * @see src/core/logging.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "assert.h"
#include "config.h"
#include "logging.h"
#include "timer.h"
#include "types_num.h"
#include "types_eig.h"

#ifdef SUPPORT_GPU
#include "../cuda/assert_cuda.h"
#include "../cuda/cuda.h"
#endif
