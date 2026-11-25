/**
 * @file section_beam.cpp
 * @brief Implements beam section reporting.
 *
 * @see src/section/section_beam.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "section_beam.h"

#include "../core/logging.h"

#include <sstream>
#include <utility>

namespace fem {
void BeamSection::info() {
    logging::info(true, "BeamSection:");
    logging::info(true, "   Material: ", (material_ ? material_->name : "-"));
    logging::info(true, "   Region  : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Profile : ", (profile_ ? profile_->name : "-"));
    logging::info(true, "   n1      : ", direction_.transpose());
}

std::string BeamSection::str() const {
    std::ostringstream os;

    os << "BeamSection: material="
       << (material_ ? material_->name : std::string("-"))
       << ", region="
       << (region_ ? region_->name : std::string("-"))
       << " ("
       << (region_ ? static_cast<int>(region_->size()) : 0)
       << ")"
       << ", profile="
       << (profile_ ? profile_->name : std::string("-"))
       << ", n1=["
       << direction_[0] << " "
       << direction_[1] << " "
       << direction_[2]
       << "]";

    return os.str();
}
// Copies the section for one compiled instance; the copy keeps the section id
Section::Ptr BeamSection::instantiate(model::ElementRegion::Ptr region,
                                      const Mat3&               rotation,
                                      const Vec3&               translation) const {
    (void) translation;
    auto copy        = std::make_shared<BeamSection>(*this);
    copy->region_    = std::move(region);
    copy->direction_ = rotation * direction_;
    return copy;
}

} // namespace fem
