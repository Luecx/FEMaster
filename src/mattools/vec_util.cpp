#include "vec_util.h"

#include "../core/core.h"

namespace fem { namespace mattools {

Vec3 normalized(Vec3 v) {
    const Precision n = v.norm();
    logging::error(n > Precision(1e-14), "Unable to normalise vector");
    return v / n;
}

Mat3 skew(const Vec3& v) {
    Mat3 S;
    S << Precision(0), -v(2),       v(1),
         v(2),        Precision(0), -v(0),
        -v(1),        v(0),        Precision(0);
    return S;
}

} } // namespace fem::mattools
