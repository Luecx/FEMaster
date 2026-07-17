#pragma once

#include "shell_material_stress.h"

namespace fem {

struct ShellMaterialStressPK2 : ShellMaterialStress {
    ShellMaterialStressPK2() = default;
    explicit ShellMaterialStressPK2(const Vec5& values);
};

} // namespace fem
