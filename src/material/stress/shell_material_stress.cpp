#include "shell_material_stress.h"

namespace fem {

ShellMaterialStress::ShellMaterialStress(const Vec5& values)
    : values_(values) {}

Precision& ShellMaterialStress::operator[](Component component) {
    return values_(static_cast<int>(component));
}

Precision ShellMaterialStress::operator[](Component component) const {
    return values_(static_cast<int>(component));
}

const Vec5& ShellMaterialStress::values() const {
    return values_;
}

Vec5& ShellMaterialStress::values() {
    return values_;
}

} // namespace fem
