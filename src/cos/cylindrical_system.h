#pragma once

#include "coordinate_system.h"
#include "../core/types_eig.h"
#include <cmath>

namespace fem {
namespace cos {

/******************************************************************************
 * @class CylindricalSystem
 * @brief Class for representing a cylindrical coordinate system.
 *
 * @details The CylindricalSystem class defines a local coordinate system based
 *          on three points in the global coordinate system: a base point, a point
 *          on the radial axis (r-axis), and a point in the theta-plane. The r-axis,
 *          theta-axis, and z-axis are orthogonalized to ensure a valid cylindrical
 *          coordinate system.
 ******************************************************************************/
class CylindricalSystem : public CoordinateSystem {
private:
    Vec3 base_point_; ///< The base point (origin) of the cylindrical coordinate system.
    Vec3 r_axis_;     ///< Unit vector defining the radial axis (r-axis).
    Vec3 theta_axis_; ///< Unit vector defining the angular direction (theta-axis).
    Vec3 z_axis_;     ///< Unit vector defining the z-axis (orthogonal to r and theta).

public:
    /**
     * @brief Constructor for the CylindricalSystem class.
     * @param base_point The origin of the cylindrical system.
     * @param r_point A point defining the r-axis (radial direction).
     * @param theta_point A point in the theta-plane.
     */
    CylindricalSystem(const std::string& name, const Vec3& base_point, const Vec3& r_point, const Vec3& theta_point)
        : CoordinateSystem(name), base_point_(base_point) {
        // Compute r-axis as the normalized direction from base_point to r_point
        r_axis_ = (r_point - base_point).normalized();

        // Compute an initial theta direction and orthogonalize it against r-axis
        Vec3 initial_theta = (theta_point - base_point).normalized();
        theta_axis_ = (initial_theta - r_axis_ * (initial_theta.dot(r_axis_))).normalized();

        // Compute z-axis as the cross product of r-axis and theta-axis
        z_axis_ = r_axis_.cross(theta_axis_).normalized();
    }

    /**
     * @brief Transform a point from global (Cartesian) to local (cylindrical) coordinates.
     * @param global_point The point in the global coordinate system (x, y, z).
     * @return The point in the cylindrical coordinate system (r, theta, z).
     */
    Vec3 to_local(const Vec3& global_point) const override {
        // Translate point to the base coordinate system
        Vec3 relative_point = global_point - base_point_;

        // Compute r as the projection onto the r-axis
        Precision r = relative_point.dot(r_axis_);

        // Compute theta as the angle in the theta-plane
        Precision theta = std::atan2(relative_point.dot(theta_axis_), relative_point.dot(r_axis_));

        // Compute z as the projection onto the z-axis
        Precision z = relative_point.dot(z_axis_);

        return Vec3(r, theta, z);
    }

    /**
     * @brief Transform a point from local (cylindrical) to global (Cartesian) coordinates.
     * @param local_point The point in the cylindrical coordinate system (r, theta, z).
     * @return The point in the global coordinate system (x, y, z).
     */
    Vec3 to_global(const Vec3& local_point) const override {
        Precision r = local_point.x();
        Precision theta = local_point.y();
        Precision z = local_point.z();

        // Convert local cylindrical coordinates to Cartesian coordinates
        Vec3 global_point = base_point_ +
                            r * r_axis_ +
                            std::cos(theta) * theta_axis_ +
                            std::sin(theta) * z_axis_ +
                            z * z_axis_;

        return global_point;
    }

    /**
     * @brief Get the cylindrical axes (unit vectors) at a specific local point.
     * @param local_point The point in the cylindrical coordinate system (r, theta, z).
     * @return A 3x3 matrix representing the axes of the cylindrical system at the given point.
     *
     * @details The returned matrix has the following columns:
     *          - Column 1: Radial direction (unit vector).
     *          - Column 2: Tangential direction (theta, unit vector).
     *          - Column 3: Z-direction (always constant, unit vector).
     */
    StaticMatrix<3, 3> get_axes(const Vec3& local_point) const override {
        Precision theta = local_point.y(); // Extract theta from the local point

        // Radial direction
        Vec3 radial = std::cos(theta) * r_axis_ + std::sin(theta) * theta_axis_;
        // Tangential direction (perpendicular to radial)
        Vec3 tangential = -std::sin(theta) * r_axis_ + std::cos(theta) * theta_axis_;
        // Z-direction (always constant)
        Vec3 z_dir = z_axis_;

        // Construct the transformation matrix
        StaticMatrix<3, 3> axes;
        axes.col(0) = radial.normalized();
        axes.col(1) = tangential.normalized();
        axes.col(2) = z_dir;

        return axes;
    }
};

} // namespace cos
} // namespace fem
