/**
 * @file pi.h
 * @brief Defines the circle constant in FEMaster's configured scalar precision.
 *
 * The math subsystem shares this compile-time constant between angle, rotation
 * and frequency calculations. The Precision alias selects the stored scalar
 * type; physical unit conversions remain responsibilities of the callers.
 *
 * @see fem::Precision
 * @author Finn Eggers
 * @date 10.11.2026
 */

#pragma once

#include "../core/types_num.h"

namespace fem::math {

// Circle constant for angles in radians and angular-frequency conversions
inline constexpr Precision pi = Precision(3.141592653589793238462643383279502884L);

} // namespace fem::math
