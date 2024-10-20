#include "quadrature.h"

#include <cmath>

namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_LINEAR or ORDER_CONSTANT
constexpr Point domain_iso_hex_points_1[] = {
    Point(0.0, 0.0, 0.0, 8.0)
};

constexpr Precision one_over_sqrt_3 = 1 / csqrt(3.0);

// ORDER_QUADRATIC or ORDER_CUBIC
constexpr Point domain_iso_hex_points_8[] = {
    Point( one_over_sqrt_3,  one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3,  one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3, -one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3,  one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3,  one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3, -one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3, 1.0)
};

// ORDER_QUARTIC or ORDER_QUINTIC
constexpr Precision coords_hex[] = { -csqrt(3.0 / 5.0), 0.0, csqrt(3.0 / 5.0) };
constexpr Precision weights_hex[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
constexpr Point domain_iso_hex_points_27[] = {
    // x = coords_hex[0]
    Point(coords_hex[0], coords_hex[0], coords_hex[0], weights_hex[0] * weights_hex[0] * weights_hex[0]),
    Point(coords_hex[0], coords_hex[0], coords_hex[1], weights_hex[0] * weights_hex[0] * weights_hex[1]),
    Point(coords_hex[0], coords_hex[0], coords_hex[2], weights_hex[0] * weights_hex[0] * weights_hex[2]),

    Point(coords_hex[0], coords_hex[1], coords_hex[0], weights_hex[0] * weights_hex[1] * weights_hex[0]),
    Point(coords_hex[0], coords_hex[1], coords_hex[1], weights_hex[0] * weights_hex[1] * weights_hex[1]),
    Point(coords_hex[0], coords_hex[1], coords_hex[2], weights_hex[0] * weights_hex[1] * weights_hex[2]),

    Point(coords_hex[0], coords_hex[2], coords_hex[0], weights_hex[0] * weights_hex[2] * weights_hex[0]),
    Point(coords_hex[0], coords_hex[2], coords_hex[1], weights_hex[0] * weights_hex[2] * weights_hex[1]),
    Point(coords_hex[0], coords_hex[2], coords_hex[2], weights_hex[0] * weights_hex[2] * weights_hex[2]),

    // x = coords_hex[1]
    Point(coords_hex[1], coords_hex[0], coords_hex[0], weights_hex[1] * weights_hex[0] * weights_hex[0]),
    Point(coords_hex[1], coords_hex[0], coords_hex[1], weights_hex[1] * weights_hex[0] * weights_hex[1]),
    Point(coords_hex[1], coords_hex[0], coords_hex[2], weights_hex[1] * weights_hex[0] * weights_hex[2]),

    Point(coords_hex[1], coords_hex[1], coords_hex[0], weights_hex[1] * weights_hex[1] * weights_hex[0]),
    Point(coords_hex[1], coords_hex[1], coords_hex[1], weights_hex[1] * weights_hex[1] * weights_hex[1]),
    Point(coords_hex[1], coords_hex[1], coords_hex[2], weights_hex[1] * weights_hex[1] * weights_hex[2]),

    Point(coords_hex[1], coords_hex[2], coords_hex[0], weights_hex[1] * weights_hex[2] * weights_hex[0]),
    Point(coords_hex[1], coords_hex[2], coords_hex[1], weights_hex[1] * weights_hex[2] * weights_hex[1]),
    Point(coords_hex[1], coords_hex[2], coords_hex[2], weights_hex[1] * weights_hex[2] * weights_hex[2]),

    // x = coords_hex[2]
    Point(coords_hex[2], coords_hex[0], coords_hex[0], weights_hex[2] * weights_hex[0] * weights_hex[0]),
    Point(coords_hex[2], coords_hex[0], coords_hex[1], weights_hex[2] * weights_hex[0] * weights_hex[1]),
    Point(coords_hex[2], coords_hex[0], coords_hex[2], weights_hex[2] * weights_hex[0] * weights_hex[2]),

    Point(coords_hex[2], coords_hex[1], coords_hex[0], weights_hex[2] * weights_hex[1] * weights_hex[0]),
    Point(coords_hex[2], coords_hex[1], coords_hex[1], weights_hex[2] * weights_hex[1] * weights_hex[1]),
    Point(coords_hex[2], coords_hex[1], coords_hex[2], weights_hex[2] * weights_hex[1] * weights_hex[2]),

    Point(coords_hex[2], coords_hex[2], coords_hex[0], weights_hex[2] * weights_hex[2] * weights_hex[0]),
    Point(coords_hex[2], coords_hex[2], coords_hex[1], weights_hex[2] * weights_hex[2] * weights_hex[1]),
    Point(coords_hex[2], coords_hex[2], coords_hex[2], weights_hex[2] * weights_hex[2] * weights_hex[2])
};


// Register the schemes
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_CONSTANT  , domain_iso_hex_points_1);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_LINEAR    , domain_iso_hex_points_1);

REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUADRATIC , domain_iso_hex_points_8);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_CUBIC     , domain_iso_hex_points_8);

REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUARTIC   , domain_iso_hex_points_27);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUINTIC   , domain_iso_hex_points_27);

} // namespace quadrature
} // namespace fem
