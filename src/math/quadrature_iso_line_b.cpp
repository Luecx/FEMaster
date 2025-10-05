/******************************************************************************
 * @file quadrature_iso_line_b.cpp
 * @brief Registers Gauss-Legendre schemes for the shifted isoparametric line.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> line_b_linear{{Point(0.5, 1.0)}};

constexpr std::array<Point, 2> line_b_quadratic{{
    Point(0.5 - 0.5 / csqrt(3.0), 0.5),
    Point(0.5 + 0.5 / csqrt(3.0), 0.5),
}};

constexpr std::array<Point, 3> line_b_quartic{{
    Point(0.5 - 0.5 * csqrt(3.0 / 5.0), 5.0 / 18.0),
    Point(0.5, 8.0 / 18.0),
    Point(0.5 + 0.5 * csqrt(3.0 / 5.0), 5.0 / 18.0),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_CONSTANT, line_b_linear);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_LINEAR, line_b_linear);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUADRATIC, line_b_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_CUBIC, line_b_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUARTIC, line_b_quartic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUINTIC, line_b_quartic);

} // namespace quadrature
} // namespace fem

