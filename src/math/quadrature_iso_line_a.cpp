#ifndef QUADRATURE_ISO_LINE_A_H
#define QUADRATURE_ISO_LINE_A_H

#include "quadrature.h"
#include "csqrt.h"

namespace fem {
namespace quadrature {

using namespace fem::math;
    
// ORDER_CONSTANT and ORDER_LINEAR
constexpr Point domain_iso_line_a_points_1[] = {
    Point(0.0, 2.0)
};

// ORDER_QUADRATIC and ORDER_CUBIC
constexpr Point domain_iso_line_a_points_2[] = {
    Point(-1.0 / csqrt(3.0), 1.0),
    Point(1.0 / csqrt(3.0), 1.0)
};

// ORDER_QUARTIC and ORDER_QUINTIC
constexpr Point domain_iso_line_a_points_3[] = {
    Point(-csqrt(3.0 / 5.0), 5.0 / 9.0),
    Point(0.0, 8.0 / 9.0),
    Point(csqrt(3.0 / 5.0), 5.0 / 9.0)
};

REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_CONSTANT, domain_iso_line_a_points_1);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_LINEAR  , domain_iso_line_a_points_1);

REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUADRATIC, domain_iso_line_a_points_2);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_CUBIC    , domain_iso_line_a_points_2);

REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUARTIC  , domain_iso_line_a_points_3);
REGISTER_SCHEME(DOMAIN_ISO_LINE_A, ORDER_QUINTIC  , domain_iso_line_a_points_3);

} // namespace quadrature
} // namespace fem

#endif // QUADRATURE_ISO_LINE_A_H