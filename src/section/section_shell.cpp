/**
 * @file section_shell.cpp
 * @brief Implements shell section matrices and reporting.
 *
 * @see src/section/section_shell.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_shell.h"

#include "../core/logging.h"

#include <sstream>

namespace fem {
StaticMatrix<6, 6> ShellSection::get_abd() {
    logging::error(material_ && material_->has_elasticity(),
                   "ShellSection requires a material with elasticity");

    return material_->elasticity()->get_abd(thickness_);
}

StaticMatrix<2, 2> ShellSection::get_shear() {
    logging::error(material_ && material_->has_elasticity(),
                   "ShellSection requires a material with elasticity");

    return material_->elasticity()->get_shear(thickness_);
}

bool ShellSection::has_density() const {
    return material_ && material_->has_density();
}

Precision ShellSection::get_density() const {
    if (has_density()) {
        return material_->get_density();
    }

    return Precision(0);
}

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material   : ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region     : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Thickness  : ", thickness_);
    logging::info(true, "   Density    : ", (has_density() ? std::to_string(get_density()) : "NO"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string ShellSection::str() const {
    std::ostringstream os;

    os << "ShellSection: t="
       << thickness_
       << ", material="
       << (material_ ? material_->name : std::string("-"))
       << ", density="
       << (has_density() ? std::to_string(get_density()) : std::string("-"))
       << ", orientation="
       << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="
       << (region_ ? region_->name : std::string("-"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")";

    return os.str();
}
} // namespace fem
