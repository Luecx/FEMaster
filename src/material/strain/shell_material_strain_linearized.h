#pragma once

#include "shell_material_strain.h"

namespace fem {

struct ShellMaterialStrainLinearized : ShellMaterialStrain {
    ShellMaterialStrainLinearized() = default;
    explicit ShellMaterialStrainLinearized(const Vec5& values);
};

} // namespace fem
