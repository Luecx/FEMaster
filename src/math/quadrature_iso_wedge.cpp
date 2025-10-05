/******************************************************************************
 * @file quadrature_iso_wedge.cpp
 * @brief Registers quadrature schemes for the isoparametric wedge.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 2> wedge_linear{{
    Point(1.0 / 3.0, 1.0 / 3.0, -csqrt(3.0) / 3.0, 1.0 / 2.0),
    Point(1.0 / 3.0, 1.0 / 3.0,  csqrt(3.0) / 3.0, 1.0 / 2.0),
}};

constexpr Precision wedge_coord = csqrt(0.6);
constexpr Precision wedge_weight_lower = 5.0 / 9.0;
constexpr Precision wedge_weight_middle = 8.0 / 9.0;

constexpr std::array<Point, 9> wedge_quadratic{{
    Point(1.0 / 6.0, 1.0 / 6.0, -wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
    Point(2.0 / 3.0, 1.0 / 6.0, -wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
    Point(1.0 / 6.0, 2.0 / 3.0, -wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
    Point(1.0 / 6.0, 1.0 / 6.0, 0.0, (1.0 / 6.0) * wedge_weight_middle),
    Point(2.0 / 3.0, 1.0 / 6.0, 0.0, (1.0 / 6.0) * wedge_weight_middle),
    Point(1.0 / 6.0, 2.0 / 3.0, 0.0, (1.0 / 6.0) * wedge_weight_middle),
    Point(1.0 / 6.0, 1.0 / 6.0, wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
    Point(2.0 / 3.0, 1.0 / 6.0, wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
    Point(1.0 / 6.0, 2.0 / 3.0, wedge_coord, (1.0 / 6.0) * wedge_weight_lower),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_CONSTANT, wedge_linear);
REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_SUPER_LINEAR, wedge_linear);
REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_SUPER_QUADRATIC, wedge_quadratic);

} // namespace quadrature
} // namespace fem

