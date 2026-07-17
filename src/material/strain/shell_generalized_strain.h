#pragma once

#include "../../core/types_eig.h"

namespace fem {

struct ShellGeneralizedStrain {
    ShellGeneralizedStrain() = default;
    explicit ShellGeneralizedStrain(const Vec8& values);

    [[nodiscard]] Vec3 membrane() const;
    [[nodiscard]] Vec3 curvature() const;
    [[nodiscard]] Vec2 transverse_shear() const;

    [[nodiscard]] const Vec8& values() const;
    [[nodiscard]] Vec8&       values();

private:
    // [epsilon_xx, epsilon_yy, gamma_xy, kappa_xx, kappa_yy, kappa_xy, gamma_xz, gamma_yz]
    Vec8 values_{Vec8::Zero()};
};

} // namespace fem
