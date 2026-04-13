#include "section_truss.h"

#include <sstream>

namespace fem {

void TrussSection::info() {
    logging::info(true, "TrussSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   A       : ", A);
}

std::string TrussSection::str() const {
    std::ostringstream os;
    os << "TrussSection: material=" << (material ? material->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", A=" << A;
    return os.str();
}

} // namespace fem
