/******************************************************************************
 * @file quadrature_iso_tri.cpp
 * @brief Registers quadrature schemes for the isoparametric triangle.
 *
 * The definitions populate the global registry via the `REGISTER_SCHEME` macro.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> tri_constant{{Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0)}};

constexpr std::array<Point, 3> tri_quadratic{{
    Point(1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
    Point(2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
    Point(1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
}};

constexpr std::array<Point, 4> tri_cubic{{
    Point(1.0 / 3.0, 1.0 / 3.0, -27.0 / 96.0),
    Point(1.0 / 5.0, 3.0 / 5.0, 25.0 / 96.0),
    Point(1.0 / 5.0, 1.0 / 5.0, 25.0 / 96.0),
    Point(3.0 / 5.0, 1.0 / 5.0, 25.0 / 96.0),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CONSTANT, tri_constant);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_LINEAR, tri_constant);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_QUADRATIC, tri_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CUBIC, tri_cubic);

} // namespace quadrature
} // namespace fem

