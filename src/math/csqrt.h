/**
 * @file csqrt.h
 * @brief Provides a constexpr square-root implementation.
 *
 * The helper is used by compile-time quadrature definitions that cannot rely on
 * `<cmath>` due to constant-evaluation restrictions.
 *
 * @see src/math/quadrature_iso_hex.cpp
 * @see src/math/quadrature_iso_wedge.cpp
 */

#pragma once

#include "../core/types_eig.h"

#include <limits>

namespace fem {
namespace math {

/**
 * @brief Performs one Newton step of the compile-time square root.
 *
 * @param x    Value whose square root is computed.
 * @param curr Current iteration estimate.
 * @param prev Previous iteration estimate.
 */
constexpr Precision csqrt_iter(Precision x, Precision curr, Precision prev) {
    return curr == prev ? curr : csqrt_iter(x, 0.5 * (curr + x / curr), curr);
}

/**
 * @brief Computes the square root using a constexpr Newton iteration.
 *
 * @param x Value whose square root is computed.
 * @return Square root of `x` or `NaN` for negative/invalid inputs.
 */
constexpr Precision csqrt(Precision x) {
    return x >= 0 && x < std::numeric_limits<Precision>::infinity()
               ? csqrt_iter(x, x, 0)
               : std::numeric_limits<Precision>::quiet_NaN();
}

} // namespace math
} // namespace fem

