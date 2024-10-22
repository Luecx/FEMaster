#ifndef QUADRATURE_ISO_LINE_B_H
#define QUADRATURE_ISO_LINE_B_H

#include "quadrature.h"
#include "csqrt.h"
namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_CONSTANT and ORDER_LINEAR
constexpr Point domain_iso_line_b_points_1[] = {
    Point(0.5, 1.0)
};

// ORDER_QUADRATIC and ORDER_CUBIC
constexpr Point domain_iso_line_b_points_2[] = {
    Point(0.5 - 0.5 / csqrt(3.0), 0.5),
    Point(0.5 + 0.5 / csqrt(3.0), 0.5)
};

// ORDER_QUARTIC and ORDER_QUINTIC
constexpr Point domain_iso_line_b_points_3[] = {
    Point(0.5 - 0.5 * csqrt(3.0 / 5.0), 5.0 / 18.0),
    Point(0.5, 8.0 / 18.0),
    Point(0.5 + 0.5 * csqrt(3.0 / 5.0), 5.0 / 18.0)
};


REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_CONSTANT, domain_iso_line_b_points_1);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_LINEAR  , domain_iso_line_b_points_1);

REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUADRATIC, domain_iso_line_b_points_2);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_CUBIC    , domain_iso_line_b_points_2);

REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUARTIC  , domain_iso_line_b_points_3);
REGISTER_SCHEME(DOMAIN_ISO_LINE_B, ORDER_QUINTIC  , domain_iso_line_b_points_3);

} // namespace quadrature
} // namespace fem

#endif // QUADRATURE_ISO_LINE_B_H
