#ifndef QUADRATURE_ISO_QUAD_H
#define QUADRATURE_ISO_QUAD_H

#include "quadrature.h"

namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_CONSTANT and ORDER_LINEAR
constexpr Point domain_iso_quad_points_1[] = {
    Point(0.0, 0.0, 0.0, 4.0)
};

// ORDER_QUADRATIC and ORDER_CUBIC
constexpr Point domain_iso_quad_points_4[] = {
    Point(csqrt(3.0) / 3.0, csqrt(3.0) / 3.0, 1.0),
    Point(-csqrt(3.0) / 3.0, csqrt(3.0) / 3.0, 1.0),
    Point(-csqrt(3.0) / 3.0, -csqrt(3.0) / 3.0, 1.0),
    Point(csqrt(3.0) / 3.0, -csqrt(3.0) / 3.0, 1.0)
};

// ORDER_QUARTIC and ORDER_QUINTIC
constexpr Point domain_iso_quad_points_9[] = {
    Point(-csqrt(0.6), -csqrt(0.6), 25.0 / 81.0),
    Point(csqrt(0.0), -csqrt(0.6), 40.0 / 81.0),
    Point(csqrt(0.6), -csqrt(0.6), 25.0 / 81.0),
    Point(-csqrt(0.6), csqrt(0.0), 40.0 / 81.0),
    Point(csqrt(0.0), csqrt(0.0), 64.0 / 81.0),
    Point(csqrt(0.6), csqrt(0.0), 40.0 / 81.0),
    Point(-csqrt(0.6), csqrt(0.6), 25.0 / 81.0),
    Point(csqrt(0.0), csqrt(0.6), 40.0 / 81.0),
    Point(csqrt(0.6), csqrt(0.6), 25.0 / 81.0)
};

REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_CONSTANT, domain_iso_quad_points_1);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_LINEAR  , domain_iso_quad_points_1);

REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUADRATIC, domain_iso_quad_points_4);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_CUBIC    , domain_iso_quad_points_4);

REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUARTIC  , domain_iso_quad_points_9);
REGISTER_SCHEME(DOMAIN_ISO_QUAD, ORDER_QUINTIC  , domain_iso_quad_points_9);

} // namespace quadrature
} // namespace fem

#endif // QUADRATURE_ISO_QUAD_H
