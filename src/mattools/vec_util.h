/**
 * @file vec_util.h
 * @brief Small vector and tensor utility helpers.
 */

#pragma once

#include "../core/types_eig.h"

namespace fem { namespace mattools {

Vec3 normalized(Vec3 v);

Mat3 skew(const Vec3& v);

} } // namespace fem::mattools
