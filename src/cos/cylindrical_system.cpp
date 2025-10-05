/******************************************************************************
 * @file cylindrical_system.cpp
 * @brief Implements the cylindrical coordinate system utilities.
 *
 * @see src/cos/cylindrical_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 ******************************************************************************/

#include "cylindrical_system.h"

namespace fem {
namespace cos {

CylindricalSystem::CylindricalSystem(const std::string& name,
                                     const Vec3& base_point,
                                     const Vec3& r_point,
                                     const Vec3& theta_point)
    : CoordinateSystem(name)
    , base_point_(base_point) {
    r_axis_ = (r_point - base_point).normalized();
    Vec3 initial_theta = (theta_point - base_point).normalized();
    theta_axis_ = (initial_theta - r_axis_ * (initial_theta.dot(r_axis_))).normalized();
    z_axis_ = r_axis_.cross(theta_axis_).normalized();
}

Vec3 CylindricalSystem::to_local(const Vec3& global_point) const {
    Vec3 relative = global_point - base_point_;
    Precision r = relative.dot(r_axis_);
    Precision theta = std::atan2(relative.dot(theta_axis_), relative.dot(r_axis_));
    Precision z = relative.dot(z_axis_);
    return Vec3(r, theta, z);
}

Vec3 CylindricalSystem::to_global(const Vec3& local_point) const {
    Precision r = local_point.x();
    Precision theta = local_point.y();
    Precision z = local_point.z();

    Vec3 radial = std::cos(theta) * r_axis_ + std::sin(theta) * theta_axis_;
    return base_point_ + r * radial + z * z_axis_;
}

Basis CylindricalSystem::get_axes(const Vec3& local_point) const {
    Precision theta = local_point.y();
    Vec3 radial = std::cos(theta) * r_axis_ + std::sin(theta) * theta_axis_;
    Vec3 tangential = -std::sin(theta) * r_axis_ + std::cos(theta) * theta_axis_;

    Basis axes;
    axes.col(0) = radial.normalized();
    axes.col(1) = tangential.normalized();
    axes.col(2) = z_axis_;
    return axes;
}

} // namespace cos
} // namespace fem
