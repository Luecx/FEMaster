/**
 * @file rectangular_system.h
 * @brief Declares a Cartesian coordinate system implementation.
 *
 * Provides utilities to construct axes from one, two, or three direction
 * vectors and exposes rotation helpers commonly used for connector
 * transformations.
 *
 * @see src/cos/coordinate_system.h
 * @see src/cos/cylindrical_system.h
 * @author Finn Eggers
 * @date 06.03.2025
 */

#pragma once

#include "coordinate_system.h"

#include <cmath>

namespace fem {
namespace cos {

/**
 * @class RectangularSystem
 * @brief Implements a right-handed rectangular coordinate system.
 */
class RectangularSystem : public CoordinateSystem {
public:
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis);
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis);
    RectangularSystem(const std::string& name, const Vec3& x_axis);

    Vec3 to_local(const Vec3& global_point) const override;
    Vec3 to_global(const Vec3& local_point) const override;
    Basis get_axes(const Vec3& local_point) const override;

    static RectangularSystem euler(Precision rot_x, Precision rot_y, Precision rot_z);

    static Basis rotation_x(Precision angle);
    static Basis rotation_y(Precision angle);
    static Basis rotation_z(Precision angle);
    static Basis rotation_x_derivative(Precision angle);
    static Basis rotation_y_derivative(Precision angle);
    static Basis rotation_z_derivative(Precision angle);

    static Basis derivative_rot_x(Precision rot_x, Precision rot_y, Precision rot_z);
    static Basis derivative_rot_y(Precision rot_x, Precision rot_y, Precision rot_z);
    static Basis derivative_rot_z(Precision rot_x, Precision rot_y, Precision rot_z);

private:
    void normalize();
    void orthogonalize();
    void compute_transformations();

    Basis local_to_global_{};
    Basis global_to_local_{};
    Vec3 x_axis_{};
    Vec3 y_axis_{};
    Vec3 z_axis_{};
};

} // namespace cos
} // namespace fem
