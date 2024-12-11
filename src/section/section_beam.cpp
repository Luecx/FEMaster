#include "section_beam.h"

namespace fem {

void BeamSection::info() {
    logging::info(true, "BeamSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Profile : ", (profile ? profile->name : "-"));
    logging::info(true, "   n1      : ", n1.transpose());
}

} // namespace fem
