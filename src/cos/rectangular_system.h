#pragma once

#include "coordinate_system.h"
#include "../core/types.h"

namespace fem {
namespace cos {

/******************************************************************************
 * @class RectangularSystem
 * @brief Class for representing a rectangular local coordinate system.
 *
 * @details The RectangularSystem class defines a local coordinate system using
 *          three orthogonal axes defined by the unit vectors of the global coordinate
 *          system. The transformation matrices `local_to_global_` and `global_to_local_`
 *          are used to map points and vectors between the local and global coordinate
 *          frames.
 ******************************************************************************/
class RectangularSystem : public CoordinateSystem {
private:
    StaticMatrix<3, 3> local_to_global_; ///< Transformation matrix from local to global coordinates.
    StaticMatrix<3, 3> global_to_local_; ///< Transformation matrix from global to local coordinates.

    Vec3 x_axis_; ///< Unit vector of the x-axis in the global coordinate system.
    Vec3 y_axis_; ///< Unit vector of the y-axis in the global coordinate system.
    Vec3 z_axis_; ///< Unit vector of the z-axis in the global coordinate system.

public:
    /**
     * @brief Constructor for the RectangularSystem class using three vectors.
     * @param x_axis Unit vector of the x-axis in the global coordinate system.
     * @param y_axis Unit vector of the y-axis in the global coordinate system.
     * @param z_axis Unit vector of the z-axis in the global coordinate system.
     */
    RectangularSystem(const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis)
        : CoordinateSystem(x_axis, y_axis, z_axis), x_axis_(x_axis), y_axis_(y_axis), z_axis_(z_axis) {
        orthogonalize();
        normalize();
        compute_transformations();
        axis_1 = x_axis_;
        axis_2 = y_axis_;
        axis_3 = z_axis_;
    }

    /**
     * @brief Constructor for the RectangularSystem class using two vectors.
     * @param x_axis Unit vector of the x-axis in the global coordinate system.
     * @param y_axis Unit vector of the y-axis in the global coordinate system.
     */
    RectangularSystem(const Vec3& x_axis, const Vec3& y_axis)
        : x_axis_(x_axis), y_axis_(y_axis), z_axis_(x_axis.cross(y_axis)) {
        orthogonalize();
        normalize();
        compute_transformations();
        axis_1 = x_axis_;
        axis_2 = y_axis_;
        axis_3 = z_axis_;
    }

    /**
     * @brief Constructor for the RectangularSystem class using one vector.
     * @param x_axis Unit vector of the x-axis in the global coordinate system.
     */
    RectangularSystem(const Vec3& x_axis)
        : x_axis_(x_axis), y_axis_(x_axis.cross(Vec3::UnitZ())), z_axis_(x_axis.cross(y_axis_)) {
        orthogonalize();
        normalize();
        compute_transformations();
        axis_1 = x_axis_;
        axis_2 = y_axis_;
        axis_3 = z_axis_;
    }

    /**
     * @brief Normalize the axes to unit vectors.
     *
     * @details This function ensures that all the axes are of unit length.
     */
    void normalize() {
        x_axis_.normalize();
        y_axis_.normalize();
        z_axis_.normalize();
    }

    /**
     * @brief Orthogonalize the coordinate system.
     *
     * @details Ensures that the z-axis is orthogonal to the x- and y-axes, and
     *          that the y-axis is adjusted to maintain orthogonality.
     */
    void orthogonalize() {
        z_axis_ = x_axis_.cross(y_axis_).normalized();
        y_axis_ = z_axis_.cross(x_axis_).normalized();
    }

    /**
     * @brief Compute the transformation matrices from the axes.
     *
     * @details Constructs the `local_to_global_` and `global_to_local_` transformation
     *          matrices using the x, y, and z axes.
     */
    void compute_transformations() {
        local_to_global_ << x_axis_, y_axis_, z_axis_;
        global_to_local_ = local_to_global_.transpose();
    }

    /**
     * @brief Transform a point from global to local coordinates.
     * @param global_point The point in the global coordinate system.
     * @return The point in the local coordinate system.
     */
    Vec3 to_local(const Vec3& global_point) const override {
        return global_to_local_ * global_point;
    }

    /**
     * @brief Transform a point from local to global coordinates.
     * @param local_point The point in the local coordinate system.
     * @return The point in the global coordinate system.
     */
    Vec3 to_global(const Vec3& local_point) const override {
        return local_to_global_ * local_point;
    }
};

} // namespace cos
} // namespace fem
