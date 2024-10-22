#include "quadrature.h"
#include <cmath>

namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_LINEAR or ORDER_CONSTANT
constexpr Point domain_iso_pyramid_points_1[] = {
    Point(0.0, 0.0, 1.0 / 4.0, 4.0 / 3.0)
};

// ORDER_QUADRATIC or ORDER_CUBIC
constexpr Precision z_lower = 0.122514822655441;
constexpr Precision z_upper = 0.544151844011225;
constexpr Precision w_lower = 0.232547451253500;
constexpr Precision w_upper = 0.100785882079825;
constexpr Precision x_upper = 0.263184055569713;
constexpr Precision x_lower = 0.506616303349787;

constexpr Point domain_iso_pyramid_points_8[] = {
    // Upper layer
    Point( x_upper,  x_upper, z_upper, w_upper),
    Point(-x_upper,  x_upper, z_upper, w_upper),
    Point(-x_upper, -x_upper, z_upper, w_upper),
    Point( x_upper, -x_upper, z_upper, w_upper),
    // Lower layer
    Point( x_lower,  x_lower, z_lower, w_lower),
    Point(-x_lower,  x_lower, z_lower, w_lower),
    Point(-x_lower, -x_lower, z_lower, w_lower),
    Point( x_lower, -x_lower, z_lower, w_lower)
};

// ORDER_QUARTIC or ORDER_QUINTIC
constexpr Precision gr_pyramid[] = {
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.228504305653967,  0.505808707853924,  0.718055741319888,
    0.228504305653967,  0.505808707853924,  0.718055741319888,
    0.228504305653967,  0.505808707853924,  0.718055741319888
};

constexpr Precision gs_pyramid[] = {
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.228504305653967,  0.505808707853924,  0.718055741319888,
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.228504305653967,  0.505808707853924,  0.718055741319888,
    -0.228504305653967, -0.505808707853924, -0.718055741319888,
    0.000000000000000,  0.000000000000000,  0.000000000000000,
    0.228504305653967,  0.505808707853924,  0.718055741319888
};

constexpr Precision gt_pyramid[] = {
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150,
    0.705002209888498,  0.347003766038352,  0.072994024073150
};

constexpr Precision weights_pyramid[] = {
    0.009244044138451, 0.045137737425885, 0.048498876871879,
    0.014790470621521, 0.072220379881415, 0.0775982029950066,
    0.009244044138451, 0.045137737425885, 0.048498876871879,
    0.014790470621521, 0.072220379881415, 0.0775982029950066,
    0.023664752994434, 0.115552607810264, 0.124157124792009,
    0.014790470621521, 0.072220379881415, 0.0775982029950066,
    0.009244044138451, 0.045137737425885, 0.048498876871879,
    0.014790470621521, 0.072220379881415, 0.0775982029950066,
    0.009244044138451, 0.045137737425885, 0.048498876871879
};

constexpr Point domain_iso_pyramid_points_27[] = {
    // Total of 27 points
    Point(gr_pyramid[0], gs_pyramid[0], gt_pyramid[0], weights_pyramid[0]),
    Point(gr_pyramid[1], gs_pyramid[1], gt_pyramid[1], weights_pyramid[1]),
    Point(gr_pyramid[2], gs_pyramid[2], gt_pyramid[2], weights_pyramid[2]),
    Point(gr_pyramid[3], gs_pyramid[3], gt_pyramid[3], weights_pyramid[3]),
    Point(gr_pyramid[4], gs_pyramid[4], gt_pyramid[4], weights_pyramid[4]),
    Point(gr_pyramid[5], gs_pyramid[5], gt_pyramid[5], weights_pyramid[5]),
    Point(gr_pyramid[6], gs_pyramid[6], gt_pyramid[6], weights_pyramid[6]),
    Point(gr_pyramid[7], gs_pyramid[7], gt_pyramid[7], weights_pyramid[7]),
    Point(gr_pyramid[8], gs_pyramid[8], gt_pyramid[8], weights_pyramid[8]),
    Point(gr_pyramid[9], gs_pyramid[9], gt_pyramid[9], weights_pyramid[9]),
    Point(gr_pyramid[10], gs_pyramid[10], gt_pyramid[10], weights_pyramid[10]),
    Point(gr_pyramid[11], gs_pyramid[11], gt_pyramid[11], weights_pyramid[11]),
    Point(gr_pyramid[12], gs_pyramid[12], gt_pyramid[12], weights_pyramid[12]),
    Point(gr_pyramid[13], gs_pyramid[13], gt_pyramid[13], weights_pyramid[13]),
    Point(gr_pyramid[14], gs_pyramid[14], gt_pyramid[14], weights_pyramid[14]),
    Point(gr_pyramid[15], gs_pyramid[15], gt_pyramid[15], weights_pyramid[15]),
    Point(gr_pyramid[16], gs_pyramid[16], gt_pyramid[16], weights_pyramid[16]),
    Point(gr_pyramid[17], gs_pyramid[17], gt_pyramid[17], weights_pyramid[17]),
    Point(gr_pyramid[18], gs_pyramid[18], gt_pyramid[18], weights_pyramid[18]),
    Point(gr_pyramid[19], gs_pyramid[19], gt_pyramid[19], weights_pyramid[19]),
    Point(gr_pyramid[20], gs_pyramid[20], gt_pyramid[20], weights_pyramid[20]),
    Point(gr_pyramid[21], gs_pyramid[21], gt_pyramid[21], weights_pyramid[21]),
    Point(gr_pyramid[22], gs_pyramid[22], gt_pyramid[22], weights_pyramid[22]),
    Point(gr_pyramid[23], gs_pyramid[23], gt_pyramid[23], weights_pyramid[23]),
    Point(gr_pyramid[24], gs_pyramid[24], gt_pyramid[24], weights_pyramid[24]),
    Point(gr_pyramid[25], gs_pyramid[25], gt_pyramid[25], weights_pyramid[25]),
    Point(gr_pyramid[26], gs_pyramid[26], gt_pyramid[26], weights_pyramid[26])
};

// Register the schemes globally
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_CONSTANT  , domain_iso_pyramid_points_1);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_LINEAR    , domain_iso_pyramid_points_1);

REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUADRATIC , domain_iso_pyramid_points_8);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_CUBIC     , domain_iso_pyramid_points_8);

REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUARTIC   , domain_iso_pyramid_points_27);
REGISTER_SCHEME(DOMAIN_ISO_PYRAMID, ORDER_QUINTIC   , domain_iso_pyramid_points_27);

} // namespace quadrature
} // namespace fem
