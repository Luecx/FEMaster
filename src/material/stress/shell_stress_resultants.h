#pragma once

#include "../../core/types_eig.h"

namespace fem {

struct ShellStressResultants {
    ShellStressResultants() = default;
    explicit ShellStressResultants(const Vec8& values);

    [[nodiscard]] Vec3 membrane() const;
    [[nodiscard]] Vec3 moments() const;
    [[nodiscard]] Vec2 transverse_shear() const;

    [[nodiscard]] const Vec8& values() const;
    [[nodiscard]] Vec8&       values();

private:
    // [Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy]
    Vec8 values_{Vec8::Zero()};
};

} // namespace fem
