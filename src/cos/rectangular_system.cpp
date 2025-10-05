/**
 * @file rectangular_system.cpp
 * @brief Implements the rectangular coordinate system utilities.
 *
 * @see src/cos/rectangular_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "rectangular_system.h"

#include <Eigen/Geometry>

namespace fem {
namespace cos {

RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis)
    : CoordinateSystem(name)
    , x_axis_(x_axis)
    , y_axis_(y_axis)
    , z_axis_(z_axis) {
    orthogonalize();
    normalize();
    compute_transformations();
}

RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis)
    : CoordinateSystem(name)
    , x_axis_(x_axis)
    , y_axis_(y_axis)
    , z_axis_(x_axis.cross(y_axis)) {
    orthogonalize();
    normalize();
    compute_transformations();
}

RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis)
    : CoordinateSystem(name), x_axis_(x_axis) {
    x_axis_.normalize();
    y_axis_ = x_axis_.unitOrthogonal();
    z_axis_ = x_axis_.cross(y_axis_).normalized();
    compute_transformations();
}


Vec3 RectangularSystem::to_local(const Vec3& global_point) const {
    return global_to_local_ * global_point;
}

Vec3 RectangularSystem::to_global(const Vec3& local_point) const {
    return local_to_global_ * local_point;
}

Basis RectangularSystem::get_axes(const Vec3& local_point) const {
    (void)local_point;
    return local_to_global_;
}

RectangularSystem RectangularSystem::euler(Precision rot_x, Precision rot_y, Precision rot_z) {
    Basis rot_mat = rotation_x(rot_x) * rotation_y(rot_y) * rotation_z(rot_z);
    return RectangularSystem("", rot_mat.col(0), rot_mat.col(1), rot_mat.col(2));
}

Basis RectangularSystem::rotation_x(Precision angle) {
    Basis rot;
    rot << 1, 0, 0,
           0, std::cos(angle), -std::sin(angle),
           0, std::sin(angle), std::cos(angle);
    return rot;
}

Basis RectangularSystem::rotation_y(Precision angle) {
    Basis rot;
    rot << std::cos(angle), 0, std::sin(angle),
           0, 1, 0,
           -std::sin(angle), 0, std::cos(angle);
    return rot;
}

Basis RectangularSystem::rotation_z(Precision angle) {
    Basis rot;
    rot << std::cos(angle), -std::sin(angle), 0,
           std::sin(angle), std::cos(angle), 0,
           0, 0, 1;
    return rot;
}

Basis RectangularSystem::rotation_x_derivative(Precision angle) {
    Basis rot;
    rot << 0, 0, 0,
           0, -std::sin(angle), -std::cos(angle),
           0, std::cos(angle), -std::sin(angle);
    return rot;
}

Basis RectangularSystem::rotation_y_derivative(Precision angle) {
    Basis rot;
    rot << -std::sin(angle), 0, std::cos(angle),
           0, 0, 0,
           -std::cos(angle), 0, -std::sin(angle);
    return rot;
}

Basis RectangularSystem::rotation_z_derivative(Precision angle) {
    Basis rot;
    rot << -std::sin(angle), -std::cos(angle), 0,
           std::cos(angle), -std::sin(angle), 0,
           0, 0, 0;
    return rot;
}

Basis RectangularSystem::derivative_rot_x(Precision rot_x, Precision rot_y, Precision rot_z) {
    return rotation_x_derivative(rot_x) * rotation_y(rot_y) * rotation_z(rot_z);
}

Basis RectangularSystem::derivative_rot_y(Precision rot_x, Precision rot_y, Precision rot_z) {
    return rotation_x(rot_x) * rotation_y_derivative(rot_y) * rotation_z(rot_z);
}

Basis RectangularSystem::derivative_rot_z(Precision rot_x, Precision rot_y, Precision rot_z) {
    return rotation_x(rot_x) * rotation_y(rot_y) * rotation_z_derivative(rot_z);
}

void RectangularSystem::normalize() {
    x_axis_.normalize();
    y_axis_.normalize();
    z_axis_.normalize();
}

void RectangularSystem::orthogonalize() {
    x_axis_.normalize();
    if (y_axis_.norm() == 0) y_axis_ = x_axis_.unitOrthogonal();
    else {
        y_axis_ = (y_axis_ - (x_axis_.dot(y_axis_) * x_axis_)).normalized();
        if (!std::isfinite(y_axis_.squaredNorm()) || y_axis_.squaredNorm() < 1e-30)
            y_axis_ = x_axis_.unitOrthogonal();
    }
    z_axis_ = x_axis_.cross(y_axis_).normalized();
}


void RectangularSystem::compute_transformations() {
    local_to_global_.col(0) = x_axis_;
    local_to_global_.col(1) = y_axis_;
    local_to_global_.col(2) = z_axis_;
    global_to_local_ = local_to_global_.transpose();
}

} // namespace cos
} // namespace fem
