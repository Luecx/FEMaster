#pragma once

#include "shell_material_strain.h"

namespace fem {

struct ShellMaterialStrainGreenLagrange : ShellMaterialStrain {
    ShellMaterialStrainGreenLagrange() = default;
    explicit ShellMaterialStrainGreenLagrange(const Vec5& values);
};

} // namespace fem
