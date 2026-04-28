/**
 * @file section_solid.cpp
 * @brief Implements solid section reporting.
 *
 * @see src/section/section_solid.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_solid.h"

#include "../core/logging.h"

#include <sstream>

namespace fem {
void SolidSection::info() {
    logging::info(true, "SolidSection:");
    logging::info(true, "   Material   : ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region     : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Orientation: ", (orientation_ ? orientation_->name : "-"));
}

std::string SolidSection::str() const {
    std::ostringstream os;

    os << "SolidSection: material="
       << (material_ ? material_->name : std::string("-"))
       << ", orientation="
       << (orientation_ ? orientation_->name : std::string("-"))
       << ", region="
       << (region_ ? region_->name : std::string("-"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")";

    return os.str();
}
} // namespace fem
