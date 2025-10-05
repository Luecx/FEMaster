/******************************************************************************
 * @file quadrature_iso_line_a.cpp
 * @brief Registers Gauss-Legendre schemes for the isoparametric line variant A.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> line_a_linear{{Point(0.0, 2.0)}};

constexpr std::array<Point, 2> line_a_quadratic{{
    Point(-1.0 / csqrt(3.0), 1.0),
    Point( 1.0 / csqrt(3.0), 1.0),
}};

constexpr std::array<Point, 3> line_a_quartic{{
    Point(-csqrt(3.0 / 5.0), 5.0 / 9.0),
    Point(0.0, 8.0 / 9.0),
    Point(csqrt(3.0 / 5.0), 5.0 / 9.0),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_CONSTANT, line_a_linear);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_LINEAR, line_a_linear);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUADRATIC, line_a_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_CUBIC, line_a_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUARTIC, line_a_quartic);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUINTIC, line_a_quartic);

} // namespace quadrature
} // namespace fem

