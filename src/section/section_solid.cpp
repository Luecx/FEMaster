#include "section_solid.h"
#include <sstream>

namespace fem {

void SolidSection::info() {
    logging::info(true, "SolidSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
}

std::string SolidSection::str() const {
    std::ostringstream os;
    os << "SolidSection: material=" << (material ? material->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")";
    return os.str();
}

} // namespace fem
