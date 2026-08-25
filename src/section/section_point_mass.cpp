/**
 * @file section_point_mass.cpp
 * @brief Implements concentrated point-mass section construction and reporting.
 *
 * @see src/section/section_point_mass.h
 * @author Finn Eggers
 * @date 25.08.2026
 */

#include "section_point_mass.h"

#include "../core/logging.h"

#include <sstream>
#include <utility>

namespace fem {

/**
 * Constructs a concentrated point-element property on an element region.
 *
 * The section is intentionally material-free. Negative values are rejected for
 * mass, inertia and stiffness because each stored coefficient contributes a
 * diagonal physical operator entry.
 *
 * @param region Element region receiving the point property.
 * @param mass Isotropic translational point mass.
 * @param rotary_inertia Rotary inertia about rx, ry and rz.
 * @param spring_constants Translational ground stiffnesses.
 * @param rotary_spring_constants Rotational ground stiffnesses.
 */
PointMassSection::PointMassSection(model::ElementRegion::Ptr region,
                                   Precision                 mass,
                                   Vec3                      rotary_inertia,
                                   Vec3                      spring_constants,
                                   Vec3                      rotary_spring_constants)
    : mass_(mass),
      rotary_inertia_(std::move(rotary_inertia)),
      spring_constants_(std::move(spring_constants)),
      rotary_spring_constants_(std::move(rotary_spring_constants)) {
    logging::error(region != nullptr,
        "PointMassSection: element region must not be null");
    logging::error(mass_ >= Precision(0),
        "PointMassSection: mass must be non-negative");
    logging::error((rotary_inertia_.array() >= Precision(0)).all(),
        "PointMassSection: rotary inertia must be non-negative");
    logging::error((spring_constants_.array() >= Precision(0)).all(),
        "PointMassSection: translational stiffness must be non-negative");
    logging::error((rotary_spring_constants_.array() >= Precision(0)).all(),
        "PointMassSection: rotational stiffness must be non-negative");

    region_ = std::move(region);
}

/**
 * Prints all concentrated point-section coefficients.
 */
void PointMassSection::info() {
    logging::info(true, "PointMassSection:");
    logging::info(true, "   Region           : ", (region_ ? region_->name : "-"));
    logging::info(true, "   Mass             : ", mass_);
    logging::info(true, "   Rotary inertia   : ", rotary_inertia_.transpose());
    logging::info(true, "   Spring           : ", spring_constants_.transpose());
    logging::info(true, "   Rotational spring: ", rotary_spring_constants_.transpose());
}

/**
 * Returns a compact one-line summary of the point property.
 *
 * @return Human-readable section description.
 */
std::string PointMassSection::str() const {
    std::ostringstream os;

    os << "PointMassSection: region="
       << (region_ ? region_->name : std::string("-"))
       << " (" << (region_ ? static_cast<int>(region_->size()) : 0) << ")"
       << ", mass=" << mass_
       << ", inertia=[" << rotary_inertia_(0) << ", " << rotary_inertia_(1) << ", " << rotary_inertia_(2) << "]"
       << ", spring=[" << spring_constants_(0) << ", " << spring_constants_(1) << ", " << spring_constants_(2) << "]"
       << ", rotspring=[" << rotary_spring_constants_(0) << ", " << rotary_spring_constants_(1) << ", " << rotary_spring_constants_(2) << "]";

    return os.str();
}

} // namespace fem
