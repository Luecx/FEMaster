/******************************************************************************
 * @file quadrature_iso_tet.cpp
 * @brief Registers quadrature schemes for the isoparametric tetrahedron.
 *
 * @see src/math/quadrature.h
 ******************************************************************************/

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> tet_linear{{Point(1.0 / 4.0, 1.0 / 4.0, 1.0 / 4.0, 1.0 / 6.0)}};

constexpr Precision tet_p1 = 1.0 / (3.0 * csqrt(5.0) - 5.0);
constexpr Precision tet_p2 = 1.0 / (csqrt(5.0) + 5.0);
constexpr std::array<Point, 4> tet_quadratic{{
    Point(tet_p2, tet_p2, tet_p2, 1.0 / 24.0),
    Point(tet_p1, tet_p2, tet_p2, 1.0 / 24.0),
    Point(tet_p2, tet_p1, tet_p2, 1.0 / 24.0),
    Point(tet_p2, tet_p2, tet_p1, 1.0 / 24.0),
}};

constexpr std::array<Point, 8> tet_cubic{{
    Point(0.0158359099, 0.3280546970, 0.3280546970, 0.138527967 / 6.0),
    Point(0.3280546970, 0.0158359099, 0.3280546970, 0.138527967 / 6.0),
    Point(0.3280546970, 0.3280546970, 0.0158359099, 0.138527967 / 6.0),
    Point(0.3280546970, 0.3280546970, 0.3280546970, 0.138527967 / 6.0),
    Point(0.6791431780, 0.1069522740, 0.1069522740, 0.111472033 / 6.0),
    Point(0.1069522740, 0.6791431780, 0.1069522740, 0.111472033 / 6.0),
    Point(0.1069522740, 0.1069522740, 0.6791431780, 0.111472033 / 6.0),
    Point(0.1069522740, 0.1069522740, 0.1069522740, 0.111472033 / 6.0),
}};

constexpr std::array<Point, 15> tet_quartic{{
    Point(0.25, 0.25, 0.25, 0.030283678097089),
    Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.006026785714286),
    Point(0.0, 1.0 / 3.0, 1.0 / 3.0, 0.006026785714286),
    Point(1.0 / 3.0, 0.0, 1.0 / 3.0, 0.006026785714286),
    Point(1.0 / 3.0, 1.0 / 3.0, 0.0, 0.006026785714286),
    Point(0.090909090909091, 0.090909090909091, 0.090909090909091, 0.011645249086029),
    Point(0.727272727272727, 0.090909090909091, 0.090909090909091, 0.011645249086029),
    Point(0.090909090909091, 0.727272727272727, 0.090909090909091, 0.011645249086029),
    Point(0.090909090909091, 0.090909090909091, 0.727272727272727, 0.011645249086029),
    Point(0.433449846426336, 0.066550153573664, 0.066550153573664, 0.010949141561386),
    Point(0.066550153573664, 0.433449846426336, 0.066550153573664, 0.010949141561386),
    Point(0.066550153573664, 0.066550153573664, 0.433449846426336, 0.010949141561386),
    Point(0.066550153573664, 0.433449846426336, 0.433449846426336, 0.010949141561386),
    Point(0.433449846426336, 0.066550153573664, 0.433449846426336, 0.010949141561386),
    Point(0.433449846426336, 0.433449846426336, 0.066550153573664, 0.010949141561386),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_CONSTANT, tet_linear);
REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_LINEAR, tet_linear);
REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_QUADRATIC, tet_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_CUBIC, tet_cubic);
REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_QUARTIC, tet_quartic);
REGISTER_SCHEME(DOMAIN_ISO_TET, ORDER_QUINTIC, tet_quartic);

} // namespace quadrature
} // namespace fem

