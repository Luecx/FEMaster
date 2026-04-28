/**
 * @file section_truss.cpp
 * @brief Implements truss section reporting.
 *
 * @see src/section/section_truss.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_truss.h"

#include "../core/logging.h"

#include <sstream>

namespace fem {
void TrussSection::info() {
    logging::info(true, "TrussSection:");
    logging::info(true, "   Material: ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region  : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Area    : ", area_);
}

std::string TrussSection::str() const {
    std::ostringstream os;

    os << "TrussSection: material="
       << (material_ ? material_->name : std::string("-"))
       << ", region="
       << (region_   ? region_->name : std::string("-"))
       << " ("
       << (region_   ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", A="
       << area_;

    return os.str();
}
} // namespace fem
