/**
 * @file cylindrical_system.h
 * @brief Declares a cylindrical coordinate system implementation.
 *
 * Converts between cylindrical and Cartesian representations using a configurable
 * base point and orientation vectors. Unlike a rectangular orientation, a
 * cylindrical system has a spatial origin; instance compilation therefore
 * rotates both axes and translates the base point.
 *
 * @see src/cos/coordinate_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "coordinate_system.h"

#include <cmath>

namespace fem {
namespace cos {

/**
 * @class CylindricalSystem
 * @brief Implements conversions between cylindrical and Cartesian frames.
 */
class CylindricalSystem : public CoordinateSystem {
public:
    CylindricalSystem(const std::string& name, const Vec3& base_point, const Vec3& r_point, const Vec3& theta_point);

    Vec3  to_local (const Vec3& global_point) const override;
    Vec3  to_global(const Vec3& local_point ) const override;
    Basis get_axes (const Vec3& local_point ) const override;

    /**
     * @brief Embeds this cylindrical frame by one rigid instance placement.
     *
     * The base point receives rotation and translation. The radial, tangential
     * and axial directions receive only the rotation, preserving all cylindrical
     * distances and angles under the rigid mapping.
     */
    Ptr transformed(const Mat3& rotation, const Vec3& translation) const override;

private:
    Vec3 base_point_{};
    Vec3 r_axis_{};
    Vec3 theta_axis_{};
    Vec3 z_axis_{};
};
} // namespace cos
} // namespace fem
