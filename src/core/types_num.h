/**
 * @file types_num.h
 * @brief Declares core numeric type aliases used throughout the project.
 *
 * The aliases centralise configuration of scalar precision and index types for
 * both CPU and GPU execution paths.
 *
 * @see src/core/types_eig.h
 * @see src/cuda/cuda_defs.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include <cstddef>
#include <cstdint>

namespace fem {

using Time  = std::uint64_t; ///< Time value expressed in milliseconds.

// Signed 64-bit type for all ids, indices, counts and sizes. Matches Eigen::Index on every platform.
using Index = std::ptrdiff_t;

static_assert(sizeof(Index) == 8, "FEMaster requires a 64-bit Index type");

#ifdef DOUBLE_PRECISION
using Precision = double; ///< CPU precision type (double).
#else
using Precision = float;  ///< CPU precision type (float).
#endif

#ifdef CUDA_DOUBLE_PRECISION
using CudaPrecision = double; ///< GPU precision type (double).
#else
using CudaPrecision = float;  ///< GPU precision type (float).
#endif
} // namespace fem
