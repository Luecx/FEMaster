#pragma once

#include "coordinate_system.h"
#include "../core/types_eig.h"

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
    Basis local_to_global_; ///< Transformation matrix from local to global coordinates.
    Basis global_to_local_; ///< Transformation matrix from global to local coordinates.

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
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis)
        : CoordinateSystem(name), x_axis_(x_axis), y_axis_(y_axis), z_axis_(z_axis) {
        orthogonalize();
        normalize();
        compute_transformations();
    }

    /**
     * @brief Constructor for the RectangularSystem class using two vectors.
     * @param x_axis Unit vector of the x-axis in the global coordinate system.
     * @param y_axis Unit vector of the y-axis in the global coordinate system.
     */
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis)
        : CoordinateSystem(name), x_axis_(x_axis), y_axis_(y_axis), z_axis_(x_axis.cross(y_axis)) {
        orthogonalize();
        normalize();
        compute_transformations();
    }

    /**
     * @brief Constructor for the RectangularSystem class using one vector.
     * @param x_axis Unit vector of the x-axis in the global coordinate system.
     */
    RectangularSystem(const std::string& name, const Vec3& x_axis)
    : CoordinateSystem(name), x_axis_(x_axis) {
        // Pick auxiliary axis: UnitZ unless x is nearly parallel to it
        Vec3 aux = (std::abs(x_axis_.normalized().dot(Vec3::UnitZ())) > 0.9)
                   ? Vec3::UnitY()
                   : Vec3::UnitZ();

        y_axis_ = x_axis_.cross(aux);
        z_axis_ = x_axis_.cross(y_axis_);

        orthogonalize();
        normalize();
        compute_transformations();
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

    Basis get_axes(const Vec3& local_point) const override {
        (void) local_point;
        return local_to_global_;
    }

    static Basis rotation_x(Precision angle) {
        Basis rot;
        rot << 1, 0, 0,
               0, std::cos(angle), -std::sin(angle),
               0, std::sin(angle), std::cos(angle);
        return rot;
    }

    static Basis rotation_y(Precision angle) {
        Basis rot;
        rot << std::cos(angle), 0, std::sin(angle),
               0, 1, 0,
               -std::sin(angle), 0, std::cos(angle);
        return rot;
    }

    static Basis rotation_z(Precision angle) {
        Basis rot;
        rot << std::cos(angle), -std::sin(angle), 0,
               std::sin(angle), std::cos(angle), 0,
               0, 0, 1;
        return rot;
    }

    static Basis rotation_x_derivative(Precision angle) {
        Basis rot;
        rot << 0, 0, 0,
               0, -std::sin(angle), -std::cos(angle),
               0,  std::cos(angle), -std::sin(angle);
        return rot;
    }

    static Basis rotation_y_derivative(Precision angle) {
        Basis rot;
        rot << -std::sin(angle), 0, std::cos(angle),
               0, 0, 0,
               -std::cos(angle), 0, -std::sin(angle);
        return rot;
    }

    static Basis rotation_z_derivative(Precision angle) {
        Basis rot;
        rot << -std::sin(angle), -std::cos(angle), 0,
               std::cos(angle), -std::sin(angle), 0,
               0, 0, 0;
        return rot;
    }

    /**
    * @brief Constructor for the RectangularSystem class using Euler-angles.
    */
    static RectangularSystem euler(Precision rot_x, Precision rot_y, Precision rot_z) {
        using std::cos;
        using std::sin;
        Basis rot_x_mat = rotation_x(rot_x);
        Basis rot_y_mat = rotation_y(rot_y);
        Basis rot_z_mat = rotation_z(rot_z);

        Basis rot_mat = rot_x_mat * rot_y_mat * rot_z_mat;

        Vec3 x_axis_ = rot_mat.col(0);
        Vec3 y_axis_ = rot_mat.col(1);
        Vec3 z_axis_ = rot_mat.col(2);
        RectangularSystem sys{"",x_axis_, y_axis_, z_axis_};
        return sys;
    }

    static Basis derivative_rot_x(Precision rot_x, Precision rot_y, Precision rot_z) {
        return rotation_x_derivative(rot_x) * rotation_y(rot_y) * rotation_z(rot_z);
    }

    static Basis derivative_rot_y(Precision rot_x, Precision rot_y, Precision rot_z) {
        return rotation_x(rot_x) * rotation_y_derivative(rot_y) * rotation_z(rot_z);
    }

    static Basis derivative_rot_z(Precision rot_x, Precision rot_y, Precision rot_z) {
        return rotation_x(rot_x) * rotation_y(rot_y) * rotation_z_derivative(rot_z);
    }
};

} // namespace cos
} // namespace fem
