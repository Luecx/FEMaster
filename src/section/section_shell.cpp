#include "section_shell.h"

namespace fem {

void ShellSection::info() {
    logging::info(true, "ShellSection:");
    logging::info(true, "   Material: ", (material ? material->name : "-"));
    logging::info(true, "   Region  : ", (region ? region->name : "-"));
    logging::info(true, "   Thickness: ", thickness);
}

} // namespace fem
