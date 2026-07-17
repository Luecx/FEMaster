#pragma once

#include "../../core/types_eig.h"

namespace fem {

struct ShellMaterialStrain {
    enum class Component : int {
        XX      = 0,
        YY      = 1,
        GammaXY = 2,
        GammaXZ = 3,
        GammaYZ = 4
    };

    ShellMaterialStrain() = default;
    explicit ShellMaterialStrain(const Vec5& values);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    [[nodiscard]] const Vec5& values() const;
    [[nodiscard]] Vec5&       values();

protected:
    Vec5 values_{Vec5::Zero()};
};

} // namespace fem
