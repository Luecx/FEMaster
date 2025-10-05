/**
 * @file quadrature_iso_quad.cpp
 * @brief Registers quadrature schemes for the isoparametric quadrilateral.
 *
 * @see src/math/quadrature.h
 */

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> quad_linear{{Point(0.0, 0.0, 0.0, 4.0)}};

constexpr std::array<Point, 4> quad_quadratic{{
    Point(fem::math::csqrt(3.0) / 3.0, fem::math::csqrt(3.0) / 3.0, 1.0),
    Point(-fem::math::csqrt(3.0) / 3.0, fem::math::csqrt(3.0) / 3.0, 1.0),
    Point(-fem::math::csqrt(3.0) / 3.0, -fem::math::csqrt(3.0) / 3.0, 1.0),
    Point(fem::math::csqrt(3.0) / 3.0, -fem::math::csqrt(3.0) / 3.0, 1.0),
}};

constexpr std::array<Point, 9> quad_quartic{{
    Point(-fem::math::csqrt(0.6), -fem::math::csqrt(0.6), 25.0 / 81.0),
    Point(0.0, -fem::math::csqrt(0.6), 40.0 / 81.0),
    Point(fem::math::csqrt(0.6), -fem::math::csqrt(0.6), 25.0 / 81.0),
    Point(-fem::math::csqrt(0.6), 0.0, 40.0 / 81.0),
    Point(0.0, 0.0, 64.0 / 81.0),
    Point(fem::math::csqrt(0.6), 0.0, 40.0 / 81.0),
    Point(-fem::math::csqrt(0.6), fem::math::csqrt(0.6), 25.0 / 81.0),
    Point(0.0, fem::math::csqrt(0.6), 40.0 / 81.0),
    Point(fem::math::csqrt(0.6), fem::math::csqrt(0.6), 25.0 / 81.0),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_CONSTANT, quad_linear);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_LINEAR, quad_linear);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUADRATIC, quad_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_CUBIC, quad_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUARTIC, quad_quartic);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUINTIC, quad_quartic);

} // namespace quadrature
} // namespace fem

