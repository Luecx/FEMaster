#include "section_solid.h.h"

namespace fem {

void SolidSection::info() {
    logging::info(true, "SolidSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
}

} // namespace fem
