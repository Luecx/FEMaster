#include "section_shell.h"
#include <sstream>

namespace fem {

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Thickness: ", thickness);
}

std::string ShellSection::str() const {
    std::ostringstream os;
    os << "ShellSection: t=" << thickness
       << ", material=" << (material ? material->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")";
    return os.str();
}

} // namespace fem
