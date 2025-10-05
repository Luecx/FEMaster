/**
 * @file quadrature_iso_pyramid.cpp
 * @brief Registers quadrature schemes for the isoparametric pyramid.
 *
 * @see src/math/quadrature.h
 */

#include "quadrature.h"

namespace fem {
namespace quadrature {

namespace {
constexpr std::array<Point, 1> pyramid_linear{{Point(0.0, 0.0, 1.0 / 4.0, 4.0 / 3.0)}};

constexpr Precision pyramid_z_lower = 0.122514822655441;
constexpr Precision pyramid_z_upper = 0.544151844011225;
constexpr Precision pyramid_weight_lower = 0.232547451253500;
constexpr Precision pyramid_weight_upper = 0.100785882079825;
constexpr Precision pyramid_x_upper = 0.263184055569713;
constexpr Precision pyramid_x_lower = 0.506616303349787;

constexpr std::array<Point, 8> pyramid_quadratic{{
    Point( pyramid_x_upper,  pyramid_x_upper, pyramid_z_upper, pyramid_weight_upper),
    Point(-pyramid_x_upper,  pyramid_x_upper, pyramid_z_upper, pyramid_weight_upper),
    Point(-pyramid_x_upper, -pyramid_x_upper, pyramid_z_upper, pyramid_weight_upper),
    Point( pyramid_x_upper, -pyramid_x_upper, pyramid_z_upper, pyramid_weight_upper),
    Point( pyramid_x_lower,  pyramid_x_lower, pyramid_z_lower, pyramid_weight_lower),
    Point(-pyramid_x_lower,  pyramid_x_lower, pyramid_z_lower, pyramid_weight_lower),
    Point(-pyramid_x_lower, -pyramid_x_lower, pyramid_z_lower, pyramid_weight_lower),
    Point( pyramid_x_lower, -pyramid_x_lower, pyramid_z_lower, pyramid_weight_lower),
}};

constexpr std::array<Point, 27> pyramid_quartic{{
    Point(-0.228504305653967, -0.228504305653967, 0.705002209888498, 0.009244044138451),
    Point(-0.505808707853924, -0.505808707853924, 0.347003766038352, 0.045137737425885),
    Point(-0.718055741319888, -0.718055741319888, 0.072994024073150, 0.048498876871879),
    Point(-0.228504305653967, 0.0, 0.705002209888498, 0.014790470621521),
    Point(-0.505808707853924, 0.0, 0.347003766038352, 0.072220379881415),
    Point(-0.718055741319888, 0.0, 0.072994024073150, 0.0775982029950066),
    Point(-0.228504305653967, 0.228504305653967, 0.705002209888498, 0.009244044138451),
    Point(-0.505808707853924, 0.505808707853924, 0.347003766038352, 0.045137737425885),
    Point(-0.718055741319888, 0.718055741319888, 0.072994024073150, 0.048498876871879),
    Point(0.0, -0.228504305653967, 0.705002209888498, 0.014790470621521),
    Point(0.0, -0.505808707853924, 0.347003766038352, 0.072220379881415),
    Point(0.0, -0.718055741319888, 0.072994024073150, 0.0775982029950066),
    Point(0.0, 0.0, 0.705002209888498, 0.023664752994434),
    Point(0.0, 0.0, 0.347003766038352, 0.115552607810264),
    Point(0.0, 0.0, 0.072994024073150, 0.124157124792009),
    Point(0.0, 0.228504305653967, 0.705002209888498, 0.014790470621521),
    Point(0.0, 0.505808707853924, 0.347003766038352, 0.072220379881415),
    Point(0.0, 0.718055741319888, 0.072994024073150, 0.0775982029950066),
    Point(0.228504305653967, -0.228504305653967, 0.705002209888498, 0.009244044138451),
    Point(0.505808707853924, -0.505808707853924, 0.347003766038352, 0.045137737425885),
    Point(0.718055741319888, -0.718055741319888, 0.072994024073150, 0.048498876871879),
    Point(0.228504305653967, 0.0, 0.705002209888498, 0.014790470621521),
    Point(0.505808707853924, 0.0, 0.347003766038352, 0.072220379881415),
    Point(0.718055741319888, 0.0, 0.072994024073150, 0.0775982029950066),
    Point(0.228504305653967, 0.228504305653967, 0.705002209888498, 0.009244044138451),
    Point(0.505808707853924, 0.505808707853924, 0.347003766038352, 0.045137737425885),
    Point(0.718055741319888, 0.718055741319888, 0.072994024073150, 0.048498876871879),
}};
} // namespace

REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_CONSTANT, pyramid_linear);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_LINEAR, pyramid_linear);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUADRATIC, pyramid_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_CUBIC, pyramid_quadratic);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUARTIC, pyramid_quartic);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUINTIC, pyramid_quartic);

} // namespace quadrature
} // namespace fem

