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
 * @brief Generalized cross-section properties used by spatial beam elements.
 *
 * A profile stores the area, centroidal bending inertia tensor, Saint-Venant
 * torsional constant and effective transverse shear areas used by B31 and B33.
 * Product of inertia follows Iyz = integral_A(y*z*dA), without a leading minus.
 * Offsets describe the shear-point and beam reference-line positions relative
 * to the section centroid.
 *
 * Named input profiles and typed geometric factories create the same object.
 * No shape geometry is retained after property evaluation. Factory methods
 * compute exact areas and bending inertias for supported sections; torsion for
 * non-circular sections uses documented engineering approximations.
 */
struct Profile : public fem::Namable {
    // Shared ownership used by model sections and named profile definitions.
    using Ptr = std::shared_ptr<Profile>;

    // Centroidal area and bending properties. Iyz uses the positive
    // integral of y*z and the section's local transverse y/z axes.
    Precision area_               = 0;
    Precision inertia_y_          = 0;
    Precision inertia_z_          = 0;
    Precision torsion_inertia_    = 0;
    Precision product_inertia_yz_ = 0;

    // Transverse shear point (SP) and beam reference line (REF) measured
    // relative to the section mid-point (SMP).
    Precision offset_y_    = 0;
    Precision offset_z_    = 0;
    Precision reference_y_ = 0;
    Precision reference_z_ = 0;

    // Effective shear areas; the constructor defaults nonpositive inputs
    // to 5A/6 for the two transverse directions.
    Precision shear_area_y_ = 0;
    Precision shear_area_z_ = 0;

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

    // Construct directly from generalized properties. The optional inertia
    // coupling, offsets and shear areas use the established FEMaster order.
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

    // Report all generalized constants and section reference offsets.
    void info();
};
} // namespace fem
