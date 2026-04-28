/**
 * @file profile.h
 * @brief Declares beam cross-section profile properties.
 *
 * A profile stores the geometric constants used by beam elements: area,
 * bending inertias, torsional inertia and optional offsets between section
 * reference axes.
 *
 * @see src/section/profile.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "../core/types_eig.h"
#include "../data/namable.h"

#include <memory>
#include <string>

namespace fem {
/**
 * @struct Profile
 * @brief Beam cross-section constants used by beam elements.
 */
struct Profile : public Namable {
    using Ptr = std::shared_ptr<Profile>; ///< Shared pointer alias for profile ownership.

    Precision area_              = 0; ///< Cross-section area.
    Precision inertia_y_         = 0; ///< Second moment of area about local y.
    Precision inertia_z_         = 0; ///< Second moment of area about local z.
    Precision torsion_inertia_   = 0; ///< Torsional inertia.
    Precision product_inertia_yz_ = 0; ///< Product of inertia: integral_A(y*z*dA).
    Precision offset_y_          = 0; ///< Shear point offset y: SP - SMP.
    Precision offset_z_          = 0; ///< Shear point offset z: SP - SMP.
    Precision reference_y_       = 0; ///< Reference point offset y: REF - SMP.
    Precision reference_z_       = 0; ///< Reference point offset z: REF - SMP.

    /**
     * @brief Creates a beam profile with geometric constants.
     *
     * @param name Profile identifier.
     * @param area Cross-section area.
     * @param inertia_y Second moment of area about local y.
     * @param inertia_z Second moment of area about local z.
     * @param torsion_inertia Torsional inertia.
     * @param product_inertia_yz Product of inertia integral_A(y*z*dA).
     * @param offset_y Shear point offset y: SP - SMP.
     * @param offset_z Shear point offset z: SP - SMP.
     * @param reference_y Reference point offset y: REF - SMP.
     * @param reference_z Reference point offset z: REF - SMP.
     */
    Profile(const std::string& name,
            Precision area,
            Precision inertia_y,
            Precision inertia_z,
            Precision torsion_inertia,
            Precision product_inertia_yz = 0,
            Precision offset_y           = 0,
            Precision offset_z           = 0,
            Precision reference_y        = 0,
            Precision reference_z        = 0);

    /**
     * @brief Prints the profile constants to the logger.
     */
    void info();
};
} // namespace fem
