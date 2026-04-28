/**
 * @file profile.cpp
 * @brief Implements beam profile storage and reporting.
 *
 * @see src/section/profile.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "profile.h"

#include "../core/logging.h"

namespace fem {
Profile::Profile(const std::string& name,
                 Precision area,
                 Precision inertia_y,
                 Precision inertia_z,
                 Precision torsion_inertia,
                 Precision product_inertia_yz,
                 Precision offset_y,
                 Precision offset_z,
                 Precision reference_y,
                 Precision reference_z)
    : Namable(name),
      area_(area),
      inertia_y_(inertia_y),
      inertia_z_(inertia_z),
      torsion_inertia_(torsion_inertia),
      product_inertia_yz_(product_inertia_yz),
      offset_y_(offset_y),
      offset_z_(offset_z),
      reference_y_(reference_y),
      reference_z_(reference_z) {}

void Profile::info() {
    logging::info(true, "Profile: ", name);
    logging::info(true, "   Area              : ", area_);
    logging::info(true, "   Inertia y         : ", inertia_y_);
    logging::info(true, "   Inertia z         : ", inertia_z_);
    logging::info(true, "   Torsion inertia   : ", torsion_inertia_);
    logging::info(true, "   Product inertia yz: ", product_inertia_yz_);
    logging::info(true, "   Offset y          : ", offset_y_);
    logging::info(true, "   Offset z          : ", offset_z_);
    logging::info(true, "   Reference y       : ", reference_y_);
    logging::info(true, "   Reference z       : ", reference_z_);
}
} // namespace fem
