#ifndef QUADRATURE_ISO_TRI_H
#define QUADRATURE_ISO_TRI_H

#include "quadrature.h"

namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_CONSTANT and ORDER_LINEAR
constexpr Point domain_iso_tri_points_1[] = {
    Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0)
};

// ORDER_QUADRATIC
constexpr Point domain_iso_tri_points_3[] = {
    Point(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
    Point(2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
    Point(1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0)
};

// ORDER_CUBIC
constexpr Point domain_iso_tri_points_4[] = {
    Point(1.0 / 3.0, 1.0 / 3.0, -27.0 / 96.0),
    Point(1.0 / 5.0, 3.0 / 5.0, 25.0 / 96.0),
    Point(1.0 / 5.0, 1.0 / 5.0, 25.0 / 96.0),
    Point(3.0 / 5.0, 1.0 / 5.0, 25.0 / 96.0)
};

REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CONSTANT, domain_iso_tri_points_1);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_LINEAR  , domain_iso_tri_points_1);

REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_QUADRATIC, domain_iso_tri_points_3);

REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CUBIC    , domain_iso_tri_points_4);

} // namespace quadrature
} // namespace fem

#endif // QUADRATURE_ISO_TRI_H
