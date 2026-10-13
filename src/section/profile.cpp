/**
 * @file profile.cpp
 * @brief Implements generalized beam properties and typed geometry factories.
 *
 * Each geometric factory evaluates area, bending inertias and a torsional
 * constant, then stores only these properties in an ordinary Profile.
 * Beam sections therefore require no geometric subtype or string-based
 * factory dispatch. Open-section torsion uses thin-wall approximations.
 *
 * @see src/section/profile.h
 * @author Finn Eggers
 * @date 28.04.2026
 */

#include "profile.h"

#include "../core/logging.h"

#include <algorithm>
#include <cmath>

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
                 Precision reference_z,
                 Precision shear_area_y,
                 Precision shear_area_z)
    : fem::Namable(name),
      area_(area),
      inertia_y_(inertia_y),
      inertia_z_(inertia_z),
      torsion_inertia_(torsion_inertia),
      product_inertia_yz_(product_inertia_yz),
      offset_y_(offset_y),
      offset_z_(offset_z),
      reference_y_(reference_y),
      reference_z_(reference_z),
      shear_area_y_(shear_area_y > Precision(0) ? shear_area_y : Precision(5) * area / Precision(6)),
      shear_area_z_(shear_area_z > Precision(0) ? shear_area_z : Precision(5) * area / Precision(6)) {
    logging::error(area_            > Precision(0),
        "PROFILE: cross-section area must be positive");
    logging::error(inertia_y_       > Precision(0),
        "PROFILE: inertia Iy must be positive");
    logging::error(inertia_z_       > Precision(0),
        "PROFILE: inertia Iz must be positive");
    logging::error(torsion_inertia_ > Precision(0),
        "PROFILE: torsional inertia must be positive");
    logging::error(shear_area_y_ > Precision(0) && shear_area_z_ > Precision(0),
        "PROFILE: effective shear areas must be positive");
    logging::error(inertia_y_ * inertia_z_ - product_inertia_yz_ * product_inertia_yz_ > Precision(0),
        "PROFILE: bending inertia matrix must be positive definite");
}


/**
 * Creates an ordinary Profile from the five directly supplied Abaqus
 * generalized-section properties. The product of inertia retains FEMaster's
 * sign convention; no geometric representation is stored.
 */
Profile::Ptr Profile::General(Precision area, Precision i11, Precision i12,
                              Precision i22, Precision torsion) {
    // Validate the finite elastic section properties before construction
    logging::error(std::isfinite(area) && area > 0,
        "PROFILE GENERAL: area must be finite and positive");
    logging::error(std::isfinite(i11),
        "PROFILE GENERAL: I11 must be finite");
    logging::error(std::isfinite(i12),
        "PROFILE GENERAL: I12 must be finite");
    logging::error(std::isfinite(i22),
        "PROFILE GENERAL: I22 must be finite");
    logging::error(std::isfinite(torsion) && torsion > 0,
        "PROFILE GENERAL: torsional constant must be finite and positive");

    return std::make_shared<Profile>("GENERAL", area, i11, i22, torsion, i12);
}

/**
 * Evaluates the exact area and centroidal bending inertias of a solid
 * rectangle and approximates its Saint-Venant torsional constant.
 *
 * The dimensions refer to the local section y/z directions.
 */
Profile::Ptr Profile::Rectangle(Precision width, Precision height) {
    // Validate physical dimensions
    logging::error(std::isfinite(width) && width > 0,
        "PROFILE RECT: width must be finite and positive");
    logging::error(std::isfinite(height) && height > 0,
        "PROFILE RECT: height must be finite and positive");

    // Compute cross-section area and centroidal inertias
    const Precision area = width * height;
    const Precision iy   = width * std::pow(height, 3) / 12;
    const Precision iz   = height * std::pow(width, 3) / 12;

    // Approximate Saint-Venant torsion of a solid rectangle
    const Precision large = std::max(width, height);
    const Precision small = std::min(width, height);
    const Precision ratio = small / large;
    const Precision j     = large * std::pow(small, 3) *
        (Precision(1) - Precision(0.63) * ratio + Precision(0.052) * std::pow(ratio, 5)) / 3;

    return std::make_shared<Profile>("RECT", area, iy, iz, j);
}

/**
 * Converts a solid circular cross-section into area, bending inertias and
 * Saint-Venant torsion. Circular symmetry makes the polar inertia equal to J.
 */
Profile::Ptr Profile::Circle(Precision radius) {
    logging::error(std::isfinite(radius) && radius > 0,
        "PROFILE CIRC: radius must be finite and positive");

    const Precision pi   = std::acos(Precision(-1));
    const Precision area = pi * radius * radius;
    const Precision iy   = pi * std::pow(radius, 4) / 4;

    return std::make_shared<Profile>("CIRC", area, iy, iy, 2 * iy);
}

/**
 * Computes the exact area, bending inertias and torsional constant of a
 * circular annulus using an outer radius and positive wall thickness.
 */
Profile::Ptr Profile::Pipe(Precision outer_radius, Precision thickness) {
    // Validate the remaining positive interior radius
    logging::error(std::isfinite(outer_radius) && outer_radius > 0,
        "PROFILE PIPE: outer radius must be finite and positive");
    logging::error(std::isfinite(thickness) && thickness > 0,
        "PROFILE PIPE: wall thickness must be finite and positive");
    logging::error(thickness < outer_radius,
        "PROFILE PIPE: thickness must be smaller than outer radius");

    const Precision pi           = std::acos(Precision(-1));
    const Precision inner_radius = outer_radius - thickness;
    const Precision area         = pi * (std::pow(outer_radius, 2) - std::pow(inner_radius, 2));
    const Precision iy           = pi * (std::pow(outer_radius, 4) - std::pow(inner_radius, 4)) / 4;

    return std::make_shared<Profile>("PIPE", area, iy, iy, 2 * iy);
}

/**
 * Creates an equal-wall box section by delegating to the four-wall form.
 * The wall convention is right, top, left and bottom.
 */
Profile::Ptr Profile::Box(Precision width, Precision height, Precision thickness) {
    return Box(width, height, thickness, thickness, thickness, thickness);
}

/**
 * Computes the centroid and exact bending inertias of a hollow rectangle
 * by summing four non-overlapping wall rectangles. The Saint-Venant torsion
 * is approximated with the thin-walled closed-cell Bredt-Batho formula.
 *
 * Wall thicknesses follow Abaqus order: right, top, left, bottom.
 */
Profile::Ptr Profile::Box(Precision width, Precision height, Precision right,
                           Precision top, Precision left, Precision bottom) {
    // Validate the exterior geometry and non-degenerate interior
    logging::error(std::isfinite(width) && width > 0,
        "PROFILE BOX: width must be finite and positive");
    logging::error(std::isfinite(height) && height > 0,
        "PROFILE BOX: height must be finite and positive");
    logging::error(std::isfinite(right) && right > 0,
        "PROFILE BOX: right wall thickness must be finite and positive");
    logging::error(std::isfinite(left) && left > 0,
        "PROFILE BOX: left wall thickness must be finite and positive");
    logging::error(std::isfinite(top) && top > 0,
        "PROFILE BOX: top wall thickness must be finite and positive");
    logging::error(std::isfinite(bottom) && bottom > 0,
        "PROFILE BOX: bottom wall thickness must be finite and positive");
    logging::error(right + left < width,
        "PROFILE BOX: side walls leave no positive interior width");
    logging::error(top + bottom < height,
        "PROFILE BOX: horizontal walls leave no positive interior height");

    // Represent the closed section as four disjoint wall rectangles.
    // Side walls span the full height; horizontal walls span the interior width.
    const Precision inner_width = width - right - left;
    const Precision a_left      = left * height;
    const Precision a_right     = right * height;
    const Precision a_top       = inner_width * top;
    const Precision a_bottom    = inner_width * bottom;
    const Precision area        = a_left + a_right + a_top + a_bottom;

    const Precision y_left      = left / 2;
    const Precision y_right     = width - right / 2;
    const Precision y_top       = left + inner_width / 2;
    const Precision y_bottom    = y_top;
    const Precision z_left      = height / 2;
    const Precision z_right     = z_left;
    const Precision z_top       = height - top / 2;
    const Precision z_bottom    = bottom / 2;

    // Determine the composite centroid using the wall areas
    const Precision cy = (a_left * y_left + a_right * y_right +
                          a_top * y_top + a_bottom * y_bottom) / area;
    const Precision cz = (a_left * z_left + a_right * z_right +
                          a_top * z_top + a_bottom * z_bottom) / area;

    // Evaluate each wall's bending inertias around the common centroid
    const Precision iy_left   = left * std::pow(height, 3) / 12 +
                                a_left * std::pow(z_left - cz, 2);
    const Precision iy_right  = right * std::pow(height, 3) / 12 +
                                a_right * std::pow(z_right - cz, 2);
    const Precision iy_top    = inner_width * std::pow(top, 3) / 12 +
                                a_top * std::pow(z_top - cz, 2);
    const Precision iy_bottom = inner_width * std::pow(bottom, 3) / 12 +
                                a_bottom * std::pow(z_bottom - cz, 2);

    const Precision iz_left   = height * std::pow(left, 3) / 12 +
                                a_left * std::pow(y_left - cy, 2);
    const Precision iz_right  = height * std::pow(right, 3) / 12 +
                                a_right * std::pow(y_right - cy, 2);
    const Precision iz_top    = top * std::pow(inner_width, 3) / 12 +
                                a_top * std::pow(y_top - cy, 2);
    const Precision iz_bottom = bottom * std::pow(inner_width, 3) / 12 +
                                a_bottom * std::pow(y_bottom - cy, 2);

    const Precision iyz_left   = a_left * (y_left - cy) * (z_left - cz);
    const Precision iyz_right  = a_right * (y_right - cy) * (z_right - cz);
    const Precision iyz_top    = a_top * (y_top - cy) * (z_top - cz);
    const Precision iyz_bottom = a_bottom * (y_bottom - cy) * (z_bottom - cz);

    // Assemble the section inertias from all four positive wall contributions
    const Precision iy  = iy_left + iy_right + iy_top + iy_bottom;
    const Precision iz  = iz_left + iz_right + iz_top + iz_bottom;
    const Precision iyz = iyz_left + iyz_right + iyz_top + iyz_bottom;

    // Estimate closed-cell torsion along the four wall midlines
    const Precision mid_width  = width - (right + left) / 2;
    const Precision mid_height = height - (top + bottom) / 2;
    const Precision j = 4 * std::pow(mid_width * mid_height, 2) /
        (mid_height / right + mid_width / top + mid_height / left + mid_width / bottom);

    return std::make_shared<Profile>("BOX", area, iy, iz, j, iyz,
                                     0, 0, width / 2 - cy, height / 2 - cz);
}

/**
 * Converts a centered I section to centroidal bending properties by summing
 * two flanges and a web. Setting one complete flange to zero produces an
 * Abaqus T section. The open-section torsion is a thin-wall estimate.
 */
Profile::Ptr Profile::I(Precision reference_height, Precision height,
                         Precision bottom_width, Precision top_width,
                         Precision bottom_thickness, Precision top_thickness,
                         Precision web_thickness) {
    // Validate section heights, flange pairs and the web
    const bool lower = bottom_width > 0 && bottom_thickness > 0;
    const bool upper = top_width > 0 && top_thickness > 0;
    logging::error(std::isfinite(reference_height),
        "PROFILE I: reference height must be finite");
    logging::error(std::isfinite(height) && height > 0,
        "PROFILE I: height must be finite and positive");
    logging::error(std::isfinite(web_thickness) && web_thickness > 0,
        "PROFILE I: web thickness must be finite and positive");
    logging::error(std::isfinite(bottom_width) && bottom_width >= 0,
        "PROFILE I: lower flange width must be finite and non-negative");
    logging::error(std::isfinite(top_width) && top_width >= 0,
        "PROFILE I: upper flange width must be finite and non-negative");
    logging::error(std::isfinite(bottom_thickness) && bottom_thickness >= 0,
        "PROFILE I: lower flange thickness must be finite and non-negative");
    logging::error(std::isfinite(top_thickness) && top_thickness >= 0,
        "PROFILE I: upper flange thickness must be finite and non-negative");
    logging::error(bottom_thickness + top_thickness < height,
        "PROFILE I: flange thicknesses leave no positive web height");
    logging::error(lower || (bottom_width == 0 && bottom_thickness == 0),
        "PROFILE I: lower flange must have both positive dimensions or both zero");
    logging::error(upper || (top_width == 0 && top_thickness == 0),
        "PROFILE I: upper flange must have both positive dimensions or both zero");
    logging::error(lower || upper,
        "PROFILE I: at least one flange is required");
    logging::error(!lower || web_thickness < bottom_width,
        "PROFILE I: web must be narrower than the lower flange");
    logging::error(!upper || web_thickness < top_width,
        "PROFILE I: web must be narrower than the upper flange");

    // Combine the disjoint bottom flange, top flange and web
    const Precision web_height = height - bottom_thickness - top_thickness;
    const Precision a1         = bottom_width * bottom_thickness;
    const Precision a2         = top_width * top_thickness;
    const Precision a3         = web_thickness * web_height;
    const Precision z1         = bottom_thickness / 2;
    const Precision z2         = height - top_thickness / 2;
    const Precision z3         = (height + bottom_thickness - top_thickness) / 2;
    const Precision area       = a1 + a2 + a3;
    const Precision cz         = (a1 * z1 + a2 * z2 + a3 * z3) / area;

    // Determine bending inertias separately for both flanges and the web.
    // Only Iy requires a parallel-axis shift because all parts are centered in y.
    const Precision iy1 = bottom_width * std::pow(bottom_thickness, 3) / 12 +
                          a1 * std::pow(z1 - cz, 2);
    const Precision iy2 = top_width * std::pow(top_thickness, 3) / 12 +
                          a2 * std::pow(z2 - cz, 2);
    const Precision iy3 = web_thickness * std::pow(web_height, 3) / 12 +
                          a3 * std::pow(z3 - cz, 2);

    const Precision iz1 = bottom_thickness * std::pow(bottom_width, 3) / 12;
    const Precision iz2 = top_thickness * std::pow(top_width, 3) / 12;
    const Precision iz3 = web_height * std::pow(web_thickness, 3) / 12;

    // Thin-wall torsion of each plate, added without warping coupling
    const Precision j1 = bottom_width * std::pow(bottom_thickness, 3) / 3;
    const Precision j2 = top_width * std::pow(top_thickness, 3) / 3;
    const Precision j3 = web_height * std::pow(web_thickness, 3) / 3;

    // The shared y-axis makes Iyz zero for centered I- and T-sections
    const Precision iy = iy1 + iy2 + iy3;
    const Precision iz = iz1 + iz2 + iz3;
    const Precision j  = j1 + j2 + j3;

    return std::make_shared<Profile>("I", area, iy, iz, j, 0,
                                     0, 0, 0, reference_height - cz);
}

/**
 * Combines two edge-aligned flanges and a web into the centroidal properties
 * of a channel section. Unequal flanges can produce a nonzero product of
 * inertia. The torsional constant uses the thin-wall open-section estimate.
 */
Profile::Ptr Profile::Channel(Precision reference_height, Precision height,
                               Precision bottom_width, Precision top_width,
                               Precision bottom_thickness, Precision top_thickness,
                               Precision web_thickness, Precision reference_width) {
    // Validate physical dimensions and reference offsets
    logging::error(std::isfinite(reference_height),
        "PROFILE CHANNEL: reference height must be finite");
    logging::error(std::isfinite(reference_width),
        "PROFILE CHANNEL: reference width must be finite");
    logging::error(std::isfinite(height) && height > 0,
        "PROFILE CHANNEL: height must be finite and positive");
    logging::error(std::isfinite(web_thickness) && web_thickness > 0,
        "PROFILE CHANNEL: web thickness must be finite and positive");
    logging::error(std::isfinite(bottom_width) && bottom_width > web_thickness,
        "PROFILE CHANNEL: lower flange must be wider than the web");
    logging::error(std::isfinite(top_width) && top_width > web_thickness,
        "PROFILE CHANNEL: upper flange must be wider than the web");
    logging::error(std::isfinite(bottom_thickness) && bottom_thickness > 0,
        "PROFILE CHANNEL: lower flange thickness must be finite and positive");
    logging::error(std::isfinite(top_thickness) && top_thickness > 0,
        "PROFILE CHANNEL: upper flange thickness must be finite and positive");
    logging::error(bottom_thickness + top_thickness < height,
        "PROFILE CHANNEL: flange thicknesses leave no positive web height");

    // Sum flange/web areas and first moments about the corner reference
    const Precision web_height = height - bottom_thickness - top_thickness;
    const Precision a1         = bottom_width * bottom_thickness;
    const Precision a2         = top_width * top_thickness;
    const Precision a3         = web_thickness * web_height;
    const Precision y1         = bottom_width / 2;
    const Precision y2         = top_width / 2;
    const Precision y3         = web_thickness / 2;
    const Precision z1         = bottom_thickness / 2;
    const Precision z2         = height - top_thickness / 2;
    const Precision z3         = (height + bottom_thickness - top_thickness) / 2;
    const Precision area       = a1 + a2 + a3;
    const Precision cy         = (a1 * y1 + a2 * y2 + a3 * y3) / area;
    const Precision cz         = (a1 * z1 + a2 * z2 + a3 * z3) / area;

    // Compute each plate's centroidal bending and coupling contributions.
    // Unlike an I-section, unequal edge-aligned flanges generally yield Iyz != 0.
    const Precision iy1 = bottom_width * std::pow(bottom_thickness, 3) / 12 +
                          a1 * std::pow(z1 - cz, 2);
    const Precision iy2 = top_width * std::pow(top_thickness, 3) / 12 +
                          a2 * std::pow(z2 - cz, 2);
    const Precision iy3 = web_thickness * std::pow(web_height, 3) / 12 +
                          a3 * std::pow(z3 - cz, 2);

    const Precision iz1 = bottom_thickness * std::pow(bottom_width, 3) / 12 +
                          a1 * std::pow(y1 - cy, 2);
    const Precision iz2 = top_thickness * std::pow(top_width, 3) / 12 +
                          a2 * std::pow(y2 - cy, 2);
    const Precision iz3 = web_height * std::pow(web_thickness, 3) / 12 +
                          a3 * std::pow(y3 - cy, 2);

    const Precision iyz1 = a1 * (y1 - cy) * (z1 - cz);
    const Precision iyz2 = a2 * (y2 - cy) * (z2 - cz);
    const Precision iyz3 = a3 * (y3 - cy) * (z3 - cz);

    // Add plate contributions to obtain the generalized section properties
    const Precision iy  = iy1 + iy2 + iy3;
    const Precision iz  = iz1 + iz2 + iz3;
    const Precision iyz = iyz1 + iyz2 + iyz3;

    const Precision j1 = bottom_width * std::pow(bottom_thickness, 3) / 3;
    const Precision j2 = top_width * std::pow(top_thickness, 3) / 3;
    const Precision j3 = web_height * std::pow(web_thickness, 3) / 3;
    const Precision j  = j1 + j2 + j3;

    return std::make_shared<Profile>("CHANNEL", area, iy, iz, j, iyz,
                                     0, 0, reference_width - cy, reference_height - cz);
}

/**
 * Combines the horizontal and vertical legs of an L section without
 * double-counting their overlap. Centroidal moments use the parallel-axis
 * theorem, and the open-section torsion remains a thin-wall estimate.
 */
Profile::Ptr Profile::Angle(Precision width, Precision height,
                             Precision bottom_thickness, Precision side_thickness) {
    // Validate finite dimensions and the non-overlapping leg decomposition
    logging::error(std::isfinite(width) && width > 0,
        "PROFILE L: width must be finite and positive");
    logging::error(std::isfinite(height) && height > 0,
        "PROFILE L: height must be finite and positive");
    logging::error(std::isfinite(bottom_thickness) && bottom_thickness > 0,
        "PROFILE L: bottom thickness must be finite and positive");
    logging::error(std::isfinite(side_thickness) && side_thickness > 0,
        "PROFILE L: side thickness must be finite and positive");
    logging::error(bottom_thickness < height,
        "PROFILE L: bottom thickness must be smaller than height");
    logging::error(side_thickness < width,
        "PROFILE L: side thickness must be smaller than width");

    // Decompose the angle into a bottom leg and the non-overlapping side leg
    const Precision side_height = height - bottom_thickness;
    const Precision a1          = width * bottom_thickness;
    const Precision a2          = side_thickness * side_height;
    const Precision y1          = width / 2;
    const Precision y2          = side_thickness / 2;
    const Precision z1          = bottom_thickness / 2;
    const Precision z2          = (height + bottom_thickness) / 2;
    const Precision area        = a1 + a2;
    const Precision cy          = (a1 * y1 + a2 * y2) / area;
    const Precision cz          = (a1 * z1 + a2 * z2) / area;

    // Evaluate bottom and side-leg bending terms around the common centroid
    const Precision iy1 = width * std::pow(bottom_thickness, 3) / 12 +
                          a1 * std::pow(z1 - cz, 2);
    const Precision iy2 = side_thickness * std::pow(side_height, 3) / 12 +
                          a2 * std::pow(z2 - cz, 2);
    const Precision iz1 = bottom_thickness * std::pow(width, 3) / 12 +
                          a1 * std::pow(y1 - cy, 2);
    const Precision iz2 = side_height * std::pow(side_thickness, 3) / 12 +
                          a2 * std::pow(y2 - cy, 2);
    const Precision iyz1 = a1 * (y1 - cy) * (z1 - cz);
    const Precision iyz2 = a2 * (y2 - cy) * (z2 - cz);

    // Assemble moments and thin-wall torsion from the two non-overlapping legs
    const Precision iy  = iy1 + iy2;
    const Precision iz  = iz1 + iz2;
    const Precision iyz = iyz1 + iyz2;
    const Precision j1  = width * std::pow(bottom_thickness, 3) / 3;
    const Precision j2  = side_height * std::pow(side_thickness, 3) / 3;
    const Precision j   = j1 + j2;

    return std::make_shared<Profile>("L", area, iy, iz, j, iyz, 0, 0, -cy, -cz);
}

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
    logging::info(true, "   Shear area y      : ", shear_area_y_);
    logging::info(true, "   Shear area z      : ", shear_area_z_);
}
} // namespace fem
