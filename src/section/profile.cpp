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
    logging::error(std::isfinite(area) && std::isfinite(i11) && std::isfinite(i12) &&
                   std::isfinite(i22) && std::isfinite(torsion),
        "PROFILE GENERAL: properties must be finite");

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
    logging::error(std::isfinite(width) && std::isfinite(height) &&
                   width > 0 && height > 0,
        "PROFILE RECT: width and height must be finite and positive");

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
    logging::error(std::isfinite(outer_radius) && std::isfinite(thickness) &&
                   outer_radius > 0 && thickness > 0 && thickness < outer_radius,
        "PROFILE PIPE: require 0 < thickness < outer radius");

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
 * by subtracting its interior from the exterior rectangle. The Saint-Venant
 * torsional constant uses the thin-walled closed-cell Bredt-Batho formula.
 *
 * Wall thicknesses follow Abaqus order: right, top, left, bottom.
 */
Profile::Ptr Profile::Box(Precision width, Precision height, Precision right,
                           Precision top, Precision left, Precision bottom) {
    // Validate the exterior geometry and non-degenerate interior
    logging::error(std::isfinite(width) && std::isfinite(height) &&
                   std::isfinite(right) && std::isfinite(top) &&
                   std::isfinite(left) && std::isfinite(bottom) &&
                   width > right + left && height > top + bottom &&
                   right > 0 && top > 0 && left > 0 && bottom > 0,
        "PROFILE BOX: positive walls must leave a positive interior");

    // Subtract the interior rectangle to obtain centroid and bending inertias
    const Precision inner_width  = width - right - left;
    const Precision inner_height = height - top - bottom;
    const Precision outer_area   = width * height;
    const Precision inner_area   = inner_width * inner_height;
    const Precision area         = outer_area - inner_area;
    const Precision hole_y       = left + inner_width / 2;
    const Precision hole_z       = bottom + inner_height / 2;
    const Precision cy           = (outer_area * width / 2 - inner_area * hole_y) / area;
    const Precision cz           = (outer_area * height / 2 - inner_area * hole_z) / area;

    const Precision iy = (width * std::pow(height, 3) -
                          inner_width * std::pow(inner_height, 3)) / 12 +
                         outer_area * std::pow(height / 2 - cz, 2) -
                         inner_area * std::pow(hole_z - cz, 2);
    const Precision iz = (height * std::pow(width, 3) -
                          inner_height * std::pow(inner_width, 3)) / 12 +
                         outer_area * std::pow(width / 2 - cy, 2) -
                         inner_area * std::pow(hole_y - cy, 2);
    const Precision iyz = outer_area * (width / 2 - cy) * (height / 2 - cz) -
                          inner_area * (hole_y - cy) * (hole_z - cz);

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
    logging::error(std::isfinite(reference_height) && std::isfinite(height) &&
                   std::isfinite(bottom_width) && std::isfinite(top_width) &&
                   std::isfinite(bottom_thickness) && std::isfinite(top_thickness) &&
                   std::isfinite(web_thickness) && height > 0 &&
                   web_thickness > 0 && bottom_width >= 0 && top_width >= 0 &&
                   bottom_thickness >= 0 && top_thickness >= 0 &&
                   bottom_thickness + top_thickness < height &&
                   (lower || (bottom_width == 0 && bottom_thickness == 0)) &&
                   (upper || (top_width == 0 && top_thickness == 0)) &&
                   (lower || upper),
        "PROFILE I: invalid web or flange dimensions");
    logging::error((!lower || web_thickness < bottom_width) &&
                   (!upper || web_thickness < top_width),
        "PROFILE I: web thickness must be smaller than each existing flange width");

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

    // The web and flanges share the centered local y-axis, hence Iyz = 0
    const Precision iy = (bottom_width * std::pow(bottom_thickness, 3) +
                          top_width * std::pow(top_thickness, 3) +
                          web_thickness * std::pow(web_height, 3)) / 12 +
                          a1 * std::pow(z1 - cz, 2) +
                          a2 * std::pow(z2 - cz, 2) +
                          a3 * std::pow(z3 - cz, 2);
    const Precision iz = (bottom_thickness * std::pow(bottom_width, 3) +
                          top_thickness * std::pow(top_width, 3) +
                          web_height * std::pow(web_thickness, 3)) / 12;
    const Precision j  = (bottom_width * std::pow(bottom_thickness, 3) +
                          top_width * std::pow(top_thickness, 3) +
                          web_height * std::pow(web_thickness, 3)) / 3;

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
    logging::error(std::isfinite(reference_height) && std::isfinite(reference_width) &&
                   std::isfinite(height) && std::isfinite(bottom_width) &&
                   std::isfinite(top_width) && std::isfinite(bottom_thickness) &&
                   std::isfinite(top_thickness) && std::isfinite(web_thickness) &&
                   height > 0 && bottom_width > web_thickness &&
                   top_width > web_thickness && web_thickness > 0 &&
                   bottom_thickness > 0 && top_thickness > 0 &&
                   bottom_thickness + top_thickness < height,
        "PROFILE CHANNEL: invalid height, flange or web dimensions");

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

    // Use the parallel-axis theorem for centroidal bending and coupling
    const Precision iy = (bottom_width * std::pow(bottom_thickness, 3) +
                          top_width * std::pow(top_thickness, 3) +
                          web_thickness * std::pow(web_height, 3)) / 12 +
                          a1 * std::pow(z1 - cz, 2) +
                          a2 * std::pow(z2 - cz, 2) +
                          a3 * std::pow(z3 - cz, 2);
    const Precision iz = (bottom_thickness * std::pow(bottom_width, 3) +
                          top_thickness * std::pow(top_width, 3) +
                          web_height * std::pow(web_thickness, 3)) / 12 +
                          a1 * std::pow(y1 - cy, 2) +
                          a2 * std::pow(y2 - cy, 2) +
                          a3 * std::pow(y3 - cy, 2);
    const Precision iyz = a1 * (y1 - cy) * (z1 - cz) +
                          a2 * (y2 - cy) * (z2 - cz) +
                          a3 * (y3 - cy) * (z3 - cz);
    const Precision j = (bottom_width * std::pow(bottom_thickness, 3) +
                         top_width * std::pow(top_thickness, 3) +
                         web_height * std::pow(web_thickness, 3)) / 3;

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
    logging::error(std::isfinite(width) && std::isfinite(height) &&
                   std::isfinite(bottom_thickness) && std::isfinite(side_thickness) &&
                   width > side_thickness && height > bottom_thickness &&
                   bottom_thickness > 0 && side_thickness > 0,
        "PROFILE L: invalid leg dimensions or thicknesses");

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

    // Recover the centroidal moments and the sign-consistent product of inertia
    const Precision iy = (width * std::pow(bottom_thickness, 3) +
                          side_thickness * std::pow(side_height, 3)) / 12 +
                          a1 * std::pow(z1 - cz, 2) + a2 * std::pow(z2 - cz, 2);
    const Precision iz = (bottom_thickness * std::pow(width, 3) +
                          side_height * std::pow(side_thickness, 3)) / 12 +
                          a1 * std::pow(y1 - cy, 2) + a2 * std::pow(y2 - cy, 2);
    const Precision iyz = a1 * (y1 - cy) * (z1 - cz) +
                          a2 * (y2 - cy) * (z2 - cz);
    const Precision j = (width * std::pow(bottom_thickness, 3) +
                         side_height * std::pow(side_thickness, 3)) / 3;

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
