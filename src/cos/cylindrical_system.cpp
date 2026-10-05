/**
 * @file cylindrical_system.cpp
 * @brief Implements cylindrical point mappings, basis evaluation and placement.
 *
 * The coordinate-system implementation orthogonalizes two defining directions,
 * uses atan2 for azimuth extraction and rotates the transverse basis by the local
 * angle. Basis columns are unit directions in global Cartesian coordinates;
 * point reconstruction also includes the stored spatial origin.
 *
 * The current to_local() implementation returns the reference-radial projection
 * rather than a radial norm. This convention is documented at its definition.
 * Rigid placement reconstructs an independent frame from its transformed origin
 * and reference directions; sections handle any material/tensor transformations.
 *
 * @see CylindricalSystem
 * @see CoordinateSystem
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "cylindrical_system.h"

#include <Eigen/Geometry>

namespace fem {
namespace cos {

/**
 * @brief Constructs a cylindrical origin and right-handed reference triad.
 *
 * Normalize r_point - base_point, remove its component from the second defining
 * direction, then complete the axial direction with a cross product. All input
 * points and resulting axes are expressed in global Cartesian coordinates.
 * The object stores its own origin and unit directions, without retaining inputs.
 *
 * Coincident points, collinear defining directions and nonfinite coordinates are
 * not checked; callers must exclude them before normalization.
 *
 * @param name Immutable identifier of the coordinate-system definition.
 * @param base_point Global Cartesian origin of the cylindrical frame.
 * @param r_point Global point defining the positive radial direction at theta = 0.
 * @param theta_point Global point whose transverse projection defines positive theta.
 */
CylindricalSystem::CylindricalSystem(const std::string& name,
                                     const Vec3& base_point,
                                     const Vec3& r_point,
                                     const Vec3& theta_point)
    : CoordinateSystem(name)
    , base_point_(base_point) {
    // Normalize the radial direction and remove it from the theta-defining direction.
    r_axis_            = (r_point       - base_point).normalized();
    Vec3 initial_theta = (theta_point   - base_point).normalized();
    theta_axis_        = (initial_theta - r_axis_ * (initial_theta.dot(r_axis_))).normalized();
    // Complete the right-handed reference triad with e_z = e_r cross e_theta.
    z_axis_ = r_axis_.cross(theta_axis_).normalized();
}

/**
 * @brief Extracts projected radial, angular and axial coordinates of a point.
 *
 * Subtract the origin and resolve the point in the reference triad. The angle is
 * atan2(q dot e_theta, q dot e_r), in radians on the atan2 branch [-pi, pi].
 * The first component is q dot e_r, not sqrt((q dot e_r)^2 + (q dot e_theta)^2).
 * Consequently this mapping is not generally the inverse of to_global(). The
 * azimuth at zero transverse distance is geometrically undefined; the code uses
 * the library's atan2 result without a special axis treatment.
 *
 * @param global_point Point in global Cartesian coordinates.
 * @return (reference-radial projection, azimuth in radians, axial projection).
 *          The coordinate-system definition is unchanged.
 */
Vec3 CylindricalSystem::to_local(const Vec3& global_point) const {
    // Resolve the origin-relative point; retain the existing radial projection convention.
    Vec3 relative   = global_point - base_point_;
    Precision r     = relative.dot(r_axis_);
    Precision theta = std::atan2(relative.dot(theta_axis_), relative.dot(r_axis_));
    Precision z     = relative.dot(z_axis_);
    return Vec3(r, theta, z);
}

/**
 * @brief Reconstructs a global point from cylindrical coordinates.
 *
 * Rotate the reference radial direction by theta and form
 * x = base_point_ + r * e_r(theta) + z * z_axis_. Negative radial values are
 * accepted algebraically. At r = 0, theta has no effect on the reconstructed point.
 * No input validation or state update is performed.
 *
 * @param local_point (r, theta, z), with theta in radians and r/z in length units.
 * @return Point in global Cartesian coordinates.
 */
Vec3 CylindricalSystem::to_global(const Vec3& local_point) const {
    // Read radial distance, azimuth and axial distance in the local representation.
    Precision r     = local_point.x();
    Precision theta = local_point.y();
    Precision z     = local_point.z();

    // Rotate the transverse direction, then restore the spatial origin and axial offset.
    Vec3 radial = std::cos(theta) * r_axis_
                + std::sin(theta) * theta_axis_;
    return base_point_ + r * radial + z * z_axis_;
}

/**
 * @brief Evaluates the cylindrical unit basis at a supplied local azimuth.
 *
 * e_r(theta) = cos(theta) * r_axis_ + sin(theta) * theta_axis_ and
 * e_theta(theta) = -sin(theta) * r_axis_ + cos(theta) * theta_axis_. The axial
 * direction is constant. Columns map local radial/tangential/axial vector
 * components to global Cartesian components; they are not coordinate-map
 * Jacobians, whose angular column would additionally contain a radial factor.
 * The caller chooses theta even on the axis; no radial singularity is evaluated.
 *
 * @param local_point Local coordinates; only the azimuth in radians is used.
 * @return Right-handed unit basis columns (e_r, e_theta, e_z) in global coordinates.
 *          The stored reference triad remains unchanged.
 */
Basis CylindricalSystem::get_axes(const Vec3& local_point) const {
    // Rotate both transverse reference directions by the local azimuth.
    Precision theta = local_point.y();
    Vec3 radial     =  std::cos(theta) * r_axis_ + std::sin(theta) * theta_axis_;
    Vec3 tangential = -std::sin(theta) * r_axis_ + std::cos(theta) * theta_axis_;

    // Store the normalized physical unit directions as global Cartesian columns.
    Basis axes;
    axes.col(0) = radial.normalized();
    axes.col(1) = tangential.normalized();
    axes.col(2) = z_axis_;
    return axes;
}

/**
 * @brief Creates an independent cylindrical definition under rigid placement.
 *
 * Apply x' = rotation * x + translation to the origin and apply only rotation to
 * the reference directions. Construct defining points one unit along each rotated
 * transverse direction so the new constructor rebuilds the same orthonormal triad.
 * The rotation must be proper orthonormal and all inputs finite; this method does
 * not validate them. Large origins may reduce the precision of defining-point
 * subtraction during reconstruction.
 *
 * @param rotation Proper orthonormal placement rotation in global coordinates.
 * @param translation Placement translation in global Cartesian coordinates.
 * @return Newly owned definition with the same name; the source is unchanged.
 */
CoordinateSystem::Ptr CylindricalSystem::transformed(const Mat3& rotation, const Vec3& translation) const {
    // Place the spatial origin and rotate the reference unit directions.
    const Vec3 base        = rotation * base_point_  + translation;
    const Vec3 radial      = rotation * r_axis_;
    const Vec3 tangential  = rotation * theta_axis_;

    // Reconstruct an independent definition from transformed global defining points.
    return std::make_shared<CylindricalSystem>(
        name,
        base,
        base + radial,
        base + tangential
    );
}
} // namespace cos
} // namespace fem
