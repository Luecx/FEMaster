/**
 * @file quadrature_iso_tri.cpp
 * @brief Registers quadrature schemes for the isoparametric triangle.
 *
 * The definitions populate the global registry via the `REGISTER_SCHEME` macro.
 *
 * @see src/math/quadrature.h
 */

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

constexpr Precision tri_quartic_a = 0.445948490915965;
constexpr Precision tri_quartic_b = 0.091576213509771;
constexpr Precision tri_quartic_c = 1.0 - 2.0 * tri_quartic_a;
constexpr Precision tri_quartic_d = 1.0 - 2.0 * tri_quartic_b;

constexpr Precision tri_quartic_wa = 0.111690794839005;
constexpr Precision tri_quartic_wb = 0.054975871827661;

constexpr std::array<Point, 6> tri_quartic{{
    Point(tri_quartic_a, tri_quartic_a, tri_quartic_wa),
    Point(tri_quartic_a, tri_quartic_c, tri_quartic_wa),
    Point(tri_quartic_c, tri_quartic_a, tri_quartic_wa),
    Point(tri_quartic_b, tri_quartic_b, tri_quartic_wb),
    Point(tri_quartic_b, tri_quartic_d, tri_quartic_wb),
    Point(tri_quartic_d, tri_quartic_b, tri_quartic_wb),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CONSTANT, tri_constant);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_LINEAR, tri_constant);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_QUADRATIC, tri_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_CUBIC, tri_cubic);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_QUARTIC, tri_quartic);
REGISTER_SCHEME(DOMAIN_ISO_TRI, ORDER_QUINTIC, tri_quartic);
} // namespace quadrature
} // namespace fem
