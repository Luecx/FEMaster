#include "shell_material_strain.h"

namespace fem {

ShellMaterialStrain::ShellMaterialStrain(const Vec5& values)
    : values_(values) {}

Precision& ShellMaterialStrain::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellMaterialStrain::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

const Vec5& ShellMaterialStrain::values() const {
    return values_;
}

Vec5& ShellMaterialStrain::values() {
    return values_;
}

} // namespace fem
