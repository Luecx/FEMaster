/******************************************************************************
 * @file quadrature_iso_hex.cpp
 * @brief Registers quadrature schemes for the isoparametric hexahedron.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> hex_linear{{Point(0.0, 0.0, 0.0, 8.0)}};

constexpr Precision one_over_sqrt_3 = csqrt(1.0 / 3.0);

constexpr std::array<Point, 8> hex_quadratic{{
    Point( one_over_sqrt_3,  one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3,  one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3, -one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point( one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3,  one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3,  one_over_sqrt_3, -one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3, -one_over_sqrt_3,  one_over_sqrt_3, 1.0),
    Point(-one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3, 1.0),
}};

constexpr std::array<Precision, 3> hex_coords{ -csqrt(3.0 / 5.0), 0.0, csqrt(3.0 / 5.0) };
constexpr std::array<Precision, 3> hex_weights{ 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

constexpr std::array<Point, 27> hex_quartic{{
    Point(hex_coords[0], hex_coords[0], hex_coords[0], hex_weights[0] * hex_weights[0] * hex_weights[0]),
    Point(hex_coords[0], hex_coords[0], hex_coords[1], hex_weights[0] * hex_weights[0] * hex_weights[1]),
    Point(hex_coords[0], hex_coords[0], hex_coords[2], hex_weights[0] * hex_weights[0] * hex_weights[2]),
    Point(hex_coords[0], hex_coords[1], hex_coords[0], hex_weights[0] * hex_weights[1] * hex_weights[0]),
    Point(hex_coords[0], hex_coords[1], hex_coords[1], hex_weights[0] * hex_weights[1] * hex_weights[1]),
    Point(hex_coords[0], hex_coords[1], hex_coords[2], hex_weights[0] * hex_weights[1] * hex_weights[2]),
    Point(hex_coords[0], hex_coords[2], hex_coords[0], hex_weights[0] * hex_weights[2] * hex_weights[0]),
    Point(hex_coords[0], hex_coords[2], hex_coords[1], hex_weights[0] * hex_weights[2] * hex_weights[1]),
    Point(hex_coords[0], hex_coords[2], hex_coords[2], hex_weights[0] * hex_weights[2] * hex_weights[2]),
    Point(hex_coords[1], hex_coords[0], hex_coords[0], hex_weights[1] * hex_weights[0] * hex_weights[0]),
    Point(hex_coords[1], hex_coords[0], hex_coords[1], hex_weights[1] * hex_weights[0] * hex_weights[1]),
    Point(hex_coords[1], hex_coords[0], hex_coords[2], hex_weights[1] * hex_weights[0] * hex_weights[2]),
    Point(hex_coords[1], hex_coords[1], hex_coords[0], hex_weights[1] * hex_weights[1] * hex_weights[0]),
    Point(hex_coords[1], hex_coords[1], hex_coords[1], hex_weights[1] * hex_weights[1] * hex_weights[1]),
    Point(hex_coords[1], hex_coords[1], hex_coords[2], hex_weights[1] * hex_weights[1] * hex_weights[2]),
    Point(hex_coords[1], hex_coords[2], hex_coords[0], hex_weights[1] * hex_weights[2] * hex_weights[0]),
    Point(hex_coords[1], hex_coords[2], hex_coords[1], hex_weights[1] * hex_weights[2] * hex_weights[1]),
    Point(hex_coords[1], hex_coords[2], hex_coords[2], hex_weights[1] * hex_weights[2] * hex_weights[2]),
    Point(hex_coords[2], hex_coords[0], hex_coords[0], hex_weights[2] * hex_weights[0] * hex_weights[0]),
    Point(hex_coords[2], hex_coords[0], hex_coords[1], hex_weights[2] * hex_weights[0] * hex_weights[1]),
    Point(hex_coords[2], hex_coords[0], hex_coords[2], hex_weights[2] * hex_weights[0] * hex_weights[2]),
    Point(hex_coords[2], hex_coords[1], hex_coords[0], hex_weights[2] * hex_weights[1] * hex_weights[0]),
    Point(hex_coords[2], hex_coords[1], hex_coords[1], hex_weights[2] * hex_weights[1] * hex_weights[1]),
    Point(hex_coords[2], hex_coords[1], hex_coords[2], hex_weights[2] * hex_weights[1] * hex_weights[2]),
    Point(hex_coords[2], hex_coords[2], hex_coords[0], hex_weights[2] * hex_weights[2] * hex_weights[0]),
    Point(hex_coords[2], hex_coords[2], hex_coords[1], hex_weights[2] * hex_weights[2] * hex_weights[1]),
    Point(hex_coords[2], hex_coords[2], hex_coords[2], hex_weights[2] * hex_weights[2] * hex_weights[2]),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_CONSTANT, hex_linear);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_LINEAR, hex_linear);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUADRATIC, hex_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_CUBIC, hex_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUARTIC, hex_quartic);
REGISTER_SCHEME(DOMAIN_ISO_HEX, ORDER_QUINTIC, hex_quartic);

} // namespace quadrature
} // namespace fem

