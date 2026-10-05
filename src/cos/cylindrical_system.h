/**
 * @file cylindrical_system.h
 * @brief Defines a cylindrical frame with a spatial origin and rotating basis.
 *
 * The coordinate-system subsystem uses CylindricalSystem to orient material axes
 * by azimuth around an axial direction. The definition stores a global origin and
 * an orthonormal reference triad; the evaluated radial and tangential directions
 * rotate with the local angle. Constitutive and tensor transformations remain
 * responsibilities of sections and other consumers of CoordinateSystem.
 *
 * Point reconstruction uses local (r, theta, z) with theta in radians. The current
 * global-to-local mapping stores a reference-radial projection as its first
 * component; it is therefore not a general inverse of point reconstruction.
 *
 * @see CylindricalSystem
 * @see CoordinateSystem
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "coordinate_system.h"

#include <cmath>

namespace fem {
namespace cos {

/**
 * @brief Owns the origin and reference triad of a cylindrical orientation.
 *
 * base_point_ is the global Cartesian origin. r_axis_ is the direction from the
 * origin to r_point. Projecting theta_point - base_point onto the plane normal to
 * r_axis_ defines theta_axis_; z_axis_ = r_axis_ cross theta_axis_ completes the
 * right-handed orthonormal triad. Inputs must be finite, the radial direction must
 * be nonzero and the two defining directions must not be collinear. Construction
 * does not validate these preconditions or provide a degeneracy fallback.
 *
 * For local angle theta in radians, get_axes() returns the radial, tangential and
 * axial unit directions as global Cartesian columns. The basis depends only on
 * theta, so the caller supplies an angle even when the physical point is on the
 * axis, where an azimuth is geometrically undefined.
 *
 * to_global() reconstructs base + r * e_r(theta) + z * e_z. to_local() returns the
 * reference-radial projection, atan2 azimuth and axial projection; its first
 * component is not the nonnegative cylindrical radius. Evaluation is read-only.
 * Rigid placement creates an independent definition with rotated directions and
 * a translated origin, preserving the name and leaving the source unchanged.
 */
class CylindricalSystem : public CoordinateSystem {
public:
    // Construct the origin and orthonormal triad from three global points.
    CylindricalSystem(const std::string& name, const Vec3& base_point, const Vec3& r_point, const Vec3& theta_point);

    // Point mappings using radians and the documented projected-radial convention.
    Vec3  to_local (const Vec3& global_point) const override;
    Vec3  to_global(const Vec3& local_point ) const override;
    // Position-dependent radial, tangential and axial columns in global coordinates.
    Basis get_axes (const Vec3& local_point ) const override;

    // Rotate the reference directions and rotate/translate the spatial origin.
    Ptr transformed(const Mat3& rotation, const Vec3& translation) const override;

private:
    // Persistent spatial origin in global Cartesian coordinates.
    Vec3 base_point_{};
    // Right-handed unit triad in global coordinates at theta = 0.
    Vec3 r_axis_{};
    Vec3 theta_axis_{};
    Vec3 z_axis_{};
};
} // namespace cos
} // namespace fem
