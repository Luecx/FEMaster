#include "section_shell.h"
#include <sstream>

namespace fem {

StaticMatrix<6, 6> ShellSection::get_abd() {
    logging::error(material && material->has_elasticity(), "ShellSection requires a material with elasticity");
    return material->elasticity()->get_abd(thickness);
}

StaticMatrix<2, 2> ShellSection::get_shear() {
    logging::error(material && material->has_elasticity(), "ShellSection requires a material with elasticity");
    return material->elasticity()->get_shear(thickness);
}

bool ShellSection::has_density() const {
    return material && material->has_density();
}

Precision ShellSection::get_density() const {
    if (material && material->has_density()) return material->get_density();
    return Precision(0);
}

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Thickness: ", thickness);
    logging::info(true, "   Density: ", (has_density() ? std::to_string(get_density()) : "NO"));
    logging::info(true, "   Orientation: ", (orientation ? orientation->name : "-"));
}

std::string ShellSection::str() const {
    std::ostringstream os;
    os << "ShellSection: t=" << thickness
       << ", material=" << (material ? material->name : std::string("-"))
       << ", density=" << (has_density() ? std::to_string(get_density()) : std::string("-"))
       << ", orientation=" << (orientation ? orientation->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")";
    return os.str();
}

} // namespace fem
