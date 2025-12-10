#include "section_point_mass.h"
#include <sstream>

namespace fem {

void PointMassSection::info() {
    logging::info(true, "PointMassSection:");
    logging::info(true, "   Material : ", (material ? material->name : "-"));
    logging::info(true, "   Region   : ", (region ? region->name : "-"));
    logging::info(true, "   Mass     : ", mass);
    logging::info(true, "   Inertia  : ", rotary_inertia.transpose());
    logging::info(true, "   Springs  : ", spring_constants.transpose());
    logging::info(true, "   Rotations: ", rotary_spring_constants.transpose());
}

std::string PointMassSection::str() const {
    std::ostringstream os;
    os << "PointMassSection: material=" << (material ? material->name : std::string("-"))
       << ", region=" << (region ? region->name : std::string("-"))
       << " (" << (region ? static_cast<int>(region->size()) : 0) << ")"
       << ", m=" << mass
       << ", J=[" << rotary_inertia[0] << " " << rotary_inertia[1] << " " << rotary_inertia[2] << "]"
       << ", k=[" << spring_constants[0] << " " << spring_constants[1] << " " << spring_constants[2] << "]"
       << ", kr=[" << rotary_spring_constants[0] << " " << rotary_spring_constants[1] << " " << rotary_spring_constants[2] << "]";
    return os.str();
}

} // namespace fem
