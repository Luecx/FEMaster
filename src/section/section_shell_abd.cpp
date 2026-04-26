#include "section_shell_abd.h"

#include <sstream>

namespace fem {

StaticMatrix<6, 6> ABDShellSection::get_abd() {
    return abd;
}

StaticMatrix<2, 2> ABDShellSection::get_shear() {
    return shear;
}

void ABDShellSection::info() {
    logging::info(true, "ABDShellSection:");
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Thickness: ", thickness);
    logging::info(true, "   Orientation: ", (orientation ? orientation->name : "-"));
}

std::string ABDShellSection::str() const {
    std::ostringstream os;
    os << "ABDShellSection: t=" << thickness
       << ", orientation=" << (orientation ? orientation->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")";
    return os.str();
}

} // namespace fem
