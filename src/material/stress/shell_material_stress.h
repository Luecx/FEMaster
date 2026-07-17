#pragma once

#include "../../core/types_eig.h"

namespace fem {

struct ShellMaterialStress {
    enum class Component : int {
        XX = 0,
        YY = 1,
        XY = 2,
        XZ = 3,
        YZ = 4
    };

    ShellMaterialStress() = default;
    explicit ShellMaterialStress(const Vec5& values);

    Precision& operator[](Component component);
    Precision  operator[](Component component) const;

    [[nodiscard]] const Vec5& values() const;
    [[nodiscard]] Vec5&       values();

protected:
    Vec5 values_{Vec5::Zero()};
};

} // namespace fem
