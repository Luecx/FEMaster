#pragma once

#include "shell_material_stress.h"

namespace fem {

struct ShellMaterialStressCauchy : ShellMaterialStress {
    ShellMaterialStressCauchy() = default;
    explicit ShellMaterialStressCauchy(const Vec5& values);
};

} // namespace fem
