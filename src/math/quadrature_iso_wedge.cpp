#include "quadrature.h"
#include <cmath>

namespace fem {
namespace quadrature {

using namespace fem::math;

// ORDER_SUPER_LINEAR
constexpr Point domain_iso_wedge_points_4[] = {
    Point(1.0 / 3.0, 1.0 / 3.0, -csqrt(3.0) / 3.0, 1.0 / 2.0),
    Point(1.0 / 3.0, 1.0 / 3.0,  csqrt(3.0) / 3.0, 1.0 / 2.0)
};

// ORDER_SUPER_QUADRATIC
constexpr Precision csqrt_0_6 = csqrt(0.6);
constexpr Precision weight_5_9 = 5.0 / 9.0;
constexpr Precision weight_8_9 = 8.0 / 9.0;

constexpr Point domain_iso_wedge_points_9[] = {
    // Lower z = -csqrt(0.6)
    Point(1.0 / 6.0, 1.0 / 6.0, -csqrt_0_6, (1.0 / 6.0) * weight_5_9),
    Point(2.0 / 3.0, 1.0 / 6.0, -csqrt_0_6, (1.0 / 6.0) * weight_5_9),
    Point(1.0 / 6.0, 2.0 / 3.0, -csqrt_0_6, (1.0 / 6.0) * weight_5_9),
    // Middle z = 0.0
    Point(1.0 / 6.0, 1.0 / 6.0, 0.0, (1.0 / 6.0) * weight_8_9),
    Point(2.0 / 3.0, 1.0 / 6.0, 0.0, (1.0 / 6.0) * weight_8_9),
    Point(1.0 / 6.0, 2.0 / 3.0, 0.0, (1.0 / 6.0) * weight_8_9),
    // Upper z = csqrt(0.6)
    Point(1.0 / 6.0, 1.0 / 6.0, csqrt_0_6, (1.0 / 6.0) * weight_5_9),
    Point(2.0 / 3.0, 1.0 / 6.0, csqrt_0_6, (1.0 / 6.0) * weight_5_9),
    Point(1.0 / 6.0, 2.0 / 3.0, csqrt_0_6, (1.0 / 6.0) * weight_5_9)
};

// Register the schemes globally
REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_CONSTANT        , domain_iso_wedge_points_4);
REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_SUPER_LINEAR    , domain_iso_wedge_points_4);
REGISTER_SCHEME(DOMAIN_ISO_WEDGE, ORDER_SUPER_QUADRATIC , domain_iso_wedge_points_9);

} // namespace quadrature
} // namespace fem
