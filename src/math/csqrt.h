//
// Created by Finn Eggers on 20.10.24.
//

#ifndef FEMASTER_CSQRT_H
#define FEMASTER_CSQRT_H

#include "../core/types_eig.h"
#include <limits>

namespace fem :: math {

Precision constexpr csqrt_iter(Precision x, Precision curr, Precision prev) {
    return curr == prev ? curr : csqrt_iter(x, 0.5 * (curr + x / curr), curr);
}
Precision constexpr csqrt(Precision x) {
    return x >= 0 && x < std::numeric_limits<Precision>::infinity()
               ? csqrt_iter(x, x, 0)
               : std::numeric_limits<Precision>::quiet_NaN();
}

}    // namespace fem::math

#endif    // FEMASTER_CSQRT_H
