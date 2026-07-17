#include "shell_generalized_strain.h"

namespace fem {

ShellGeneralizedStrain::ShellGeneralizedStrain(const Vec8& values)
    : values_(values) {}

Vec3 ShellGeneralizedStrain::membrane() const {
    return values_.template segment<3>(0);
}

Vec3 ShellGeneralizedStrain::curvature() const {
    return values_.template segment<3>(3);
}

Vec2 ShellGeneralizedStrain::transverse_shear() const {
    return values_.template segment<2>(6);
}

const Vec8& ShellGeneralizedStrain::values() const {
    return values_;
}

Vec8& ShellGeneralizedStrain::values() {
    return values_;
}

} // namespace fem
