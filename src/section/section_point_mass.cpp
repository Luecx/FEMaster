#include "section_point_mass.h"

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

} // namespace fem
