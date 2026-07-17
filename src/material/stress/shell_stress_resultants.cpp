#include "shell_stress_resultants.h"

namespace fem {

ShellStressResultants::ShellStressResultants(const Vec8& values)
    : values_(values) {}

Vec3 ShellStressResultants::membrane() const {
    return values_.template segment<3>(0);
}

Vec3 ShellStressResultants::moments() const {
    return values_.template segment<3>(3);
}

Vec2 ShellStressResultants::transverse_shear() const {
    return values_.template segment<2>(6);
}

const Vec8& ShellStressResultants::values() const {
    return values_;
}

Vec8& ShellStressResultants::values() {
    return values_;
}

} // namespace fem
