#include "section_beam.h"
#include <sstream>

namespace fem {

void BeamSection::info() {
    logging::info(true, "BeamSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Profile : ", (profile ? profile->name : "-"));
    logging::info(true, "   n1      : ", n1.transpose());
}

std::string BeamSection::str() const {
    std::ostringstream os;
    os << "BeamSection: material=" << (material ? material->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", profile=" << (profile ? profile->name : std::string("-"))
       << ", n1=[" << n1[0] << " " << n1[1] << " " << n1[2] << "]";
    return os.str();
}

} // namespace fem
