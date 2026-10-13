/**
 * @file profile.h
 * @brief Declares generalized beam properties and typed geometry factories.
 *
 * A Profile stores only the geometric constants used by beam elements:
 * area, bending inertias, Saint-Venant torsion, shear areas and offsets.
 * Named profiles and geometrically generated profiles have identical storage.
 * Static factory methods compute these properties from standard cross-section
 * dimensions without introducing separate geometric profile classes.
 *
 * @see src/section/profile.cpp
 * @author Finn Eggers
 * @date 28.04.2026
 */

#pragma once

#include "../core/types_eig.h"
#include "../core/namable.h"

#include <memory>
#include <string>

namespace fem {
/**
 * @struct Profile
 * @brief Beam cross-section constants used by beam elements.
 */
struct Profile : public fem::Namable {
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
    Precision shear_area_y_      = 0; ///< Effective transverse shear area in local y.
    Precision shear_area_z_      = 0; ///< Effective transverse shear area in local z.

    // Construct generalized section constants from explicit properties or
    // standard cross-section dimensions. Factories return ordinary Profiles;
    // no geometric subtype is stored after evaluating the section properties.
    static Ptr General  (Precision area, Precision i11, Precision i12, Precision i22, Precision torsion);
    static Ptr Rectangle(Precision width, Precision height);
    static Ptr Circle   (Precision radius);
    static Ptr Pipe     (Precision outer_radius, Precision thickness);

    // Box with equal or individual wall thicknesses (right, top, left, bottom)
    static Ptr Box(Precision width, Precision height, Precision thickness);
    static Ptr Box(Precision width, Precision height, Precision right, Precision top,
                   Precision left, Precision bottom);

    // Abaqus I/CHANNEL dimensions: reference height, total height, flange
    // widths, flange thicknesses and web thickness. An I profile may omit one
    // complete flange (width and thickness both zero) to describe a T profile.
    static Ptr I(Precision reference_height, Precision height, Precision bottom_width,
                 Precision top_width, Precision bottom_thickness, Precision top_thickness,
                 Precision web_thickness);
    static Ptr Channel(Precision reference_height, Precision height, Precision bottom_width,
                       Precision top_width, Precision bottom_thickness, Precision top_thickness,
                       Precision web_thickness, Precision reference_width);
    static Ptr Angle(Precision width, Precision height, Precision bottom_thickness,
                     Precision side_thickness);

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
     * @param shear_area_y Effective transverse shear area in local y; non-positive selects 5A/6.
     * @param shear_area_z Effective transverse shear area in local z; non-positive selects 5A/6.
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
            Precision reference_z        = 0,
            Precision shear_area_y       = 0,
            Precision shear_area_z       = 0);

    /**
     * @brief Prints the profile constants to the logger.
     */
    void info();
};
} // namespace fem
