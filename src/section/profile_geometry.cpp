#include "profile_geometry.h"

#include "../core/logging.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <vector>

namespace fem {
namespace {

struct Rectangle {
    Precision width, height, y, z;
};

struct Properties {
    Precision area = 0;
    Precision iy = 0;
    Precision iz = 0;
    Precision iyz = 0;
    Precision cy = 0;
    Precision cz = 0;
};

// Exact area and centroidal second moments for disjoint, axis-aligned
// rectangles. A negative rectangle subtracts a hole.
Properties rectangles(std::initializer_list<Rectangle> parts) {
    Properties p{};
    for (const auto& r : parts) {
        const Precision a = r.width * r.height;
        p.area += a;
        p.cy += a * r.y;
        p.cz += a * r.z;
    }
    logging::error(p.area > 0, "BEAM SECTION: cross-section area must be positive");
    p.cy /= p.area;
    p.cz /= p.area;

    for (const auto& r : parts) {
        const Precision a = r.width * r.height;
        const Precision dy = r.y - p.cy;
        const Precision dz = r.z - p.cz;
        p.iy += r.width * r.height * r.height * r.height / 12 + a * dz * dz;
        p.iz += r.height * r.width * r.width * r.width / 12 + a * dy * dy;
        p.iyz += a * dy * dz;
    }
    return p;
}

Precision rectangle_torsion(Precision a, Precision b) {
    const Precision width = std::max(a, b);
    const Precision depth = std::min(a, b);
    const Precision ratio = depth / width;
    // Saint-Venant torsion approximation for a solid rectangle.
    return width * depth * depth * depth *
        (Precision(1) - Precision(0.63) * ratio +
         Precision(0.052) * std::pow(ratio, 5)) / Precision(3);
}

void validate_dimensions(const std::string& shape, const std::vector<Precision>& d) {
    for (std::size_t i = 0; i < d.size(); ++i) {
        // I/CHANNEL section origin offsets l and o may legitimately be zero.
        const bool offset = (shape == "I" || shape == "CHANNEL") &&
                            (i == 0 || (shape == "CHANNEL" && i == 7));
        logging::error(std::isfinite(d[i]) && (offset || d[i] > 0),
                       "BEAM SECTION: dimensions must be finite, lengths/thicknesses positive");
    }
}

} // namespace

Profile::Ptr make_beam_geometry_profile(const std::string& name,
                                        const std::string& shape,
                                        const std::vector<Precision>& d) {
    validate_dimensions(shape, d);

    const Precision pi = std::acos(Precision(-1));

    if (shape == "RECT") {
        logging::error(d.size() == 2, "BEAM SECTION RECT: expected a, b");
        const Precision a = d[0], b = d[1];
        return std::make_shared<Profile>(
            name, a * b, a * b * b * b / 12, b * a * a * a / 12,
            rectangle_torsion(a, b));
    }
    if (shape == "CIRC") {
        logging::error(d.size() == 1, "BEAM SECTION CIRC: expected radius");
        const Precision r = d[0];
        const Precision iy = pi * std::pow(r, 4) / 4;
        return std::make_shared<Profile>(name, pi * r * r, iy, iy, 2 * iy);
    }
    if (shape == "PIPE") {
        logging::error(d.size() == 2, "BEAM SECTION PIPE: expected outer radius, thickness");
        const Precision r = d[0], t = d[1];
        logging::error(t < r, "BEAM SECTION PIPE: thickness must be smaller than outer radius");
        const Precision inner = r - t;
        const Precision iy = pi * (std::pow(r, 4) - std::pow(inner, 4)) / 4;
        return std::make_shared<Profile>(
            name, pi * (r * r - inner * inner), iy, iy, 2 * iy);
    }
    if (shape == "BOX") {
        logging::error(d.size() == 3 || d.size() == 6,
                       "BEAM SECTION BOX: expected a, b, t or a, b, t1, t2, t3, t4");
        const Precision a = d[0], b = d[1];
        const Precision t1 = d[2], t2 = d.size() == 3 ? d[2] : d[3];
        const Precision t3 = d.size() == 3 ? d[2] : d[4];
        const Precision t4 = d.size() == 3 ? d[2] : d[5];
        // t1=bottom, t2=right, t3=top, t4=left (Abaqus wall order).
        const Precision inner_y = a - t2 - t4;
        const Precision inner_z = b - t1 - t3;
        logging::error(inner_y > 0 && inner_z > 0,
                       "BEAM SECTION BOX: wall thicknesses leave no interior");
        const Properties p = rectangles({
            {a, b, a / 2, b / 2},
            {-inner_y, inner_z, t4 + inner_y / 2, t1 + inner_z / 2}
        });
        const Precision mid_y = a - (t2 + t4) / 2;
        const Precision mid_z = b - (t1 + t3) / 2;
        // Thin-walled Bredt-Batho closed-cell torsion estimate.
        const Precision j = 4 * std::pow(mid_y * mid_z, 2) /
            (mid_y / t1 + mid_z / t2 + mid_y / t3 + mid_z / t4);
        return std::make_shared<Profile>(
            name, p.area, p.iy, p.iz, j, p.iyz,
            0, 0, a / 2 - p.cy, b / 2 - p.cz);
    }
    if (shape == "I" || shape == "CHANNEL") {
        logging::error(d.size() == (shape == "I" ? 7u : 8u),
                       "BEAM SECTION I/CHANNEL: expected l, h, b1, b2, t1, t2, t3 [, o]");
        const Precision l = d[0], h = d[1], b1 = d[2], b2 = d[3];
        const Precision t1 = d[4], t2 = d[5], t3 = d[6];
        logging::error(t1 + t2 < h, "BEAM SECTION I/CHANNEL: flange thicknesses exceed height");
        logging::error(t3 < std::min(b1, b2), "BEAM SECTION I/CHANNEL: web too thick");
        const Precision lower_y = shape == "I" ? 0 : b1 / 2;
        const Precision upper_y = shape == "I" ? 0 : b2 / 2;
        const Precision web_y = shape == "I" ? 0 : t3 / 2;
        const Properties p = rectangles({
            {b1, t1, lower_y, t1 / 2},
            {b2, t2, upper_y, h - t2 / 2},
            {t3, h - t1 - t2, web_y, (h + t1 - t2) / 2}
        });
        // Open-section thin-wall Saint-Venant estimate. Warping and exact
        // shear-center corrections are not represented by Profile.
        const Precision j = (b1 * std::pow(t1, 3) +
                             b2 * std::pow(t2, 3) +
                             (h - t1 - t2) * std::pow(t3, 3)) / 3;
        const Precision ref_y = (shape == "CHANNEL" ? d[7] : Precision(0)) - p.cy;
        return std::make_shared<Profile>(
            name, p.area, p.iy, p.iz, j, p.iyz, 0, 0, ref_y, l - p.cz);
    }
    if (shape == "L") {
        logging::error(d.size() == 4, "BEAM SECTION L: expected a, b, t1, t2");
        const Precision a = d[0], b = d[1], t1 = d[2], t2 = d[3];
        logging::error(t1 < b && t2 < a, "BEAM SECTION L: thicknesses exceed legs");
        const Properties p = rectangles({
            {a, t1, a / 2, t1 / 2},
            {t2, b - t1, t2 / 2, (b + t1) / 2}
        });
        const Precision j = (a * std::pow(t1, 3) + (b - t1) * std::pow(t2, 3)) / 3;
        return std::make_shared<Profile>(name, p.area, p.iy, p.iz, j, p.iyz,
                                         0, 0, -p.cy, -p.cz);
    }

    logging::error(false, "BEAM SECTION: unsupported geometric SECTION=", shape);
    return nullptr;
}

} // namespace fem
