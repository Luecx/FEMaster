/**
 * @file rectangular_system.h
 * @brief Defines Cartesian orientations and analytical rotation derivatives.
 *
 * The coordinate-system subsystem uses RectangularSystem for a constant local
 * basis without a spatial origin. Constructors complete a right-handed triad
 * from one or more supplied directions and cache both component transformations.
 * Consumers remain responsible for transforming their own tensors and DOFs.
 *
 * Rotation helpers use radians and compose R = Rx(rot_x) * Ry(rot_y) * Rz(rot_z).
 * Analytical derivatives differentiate this same product for orientation-sensitive
 * finite-element calculations. Rigid instance placement rotates the definition;
 * translation does not affect an orientation without a spatial origin.
 *
 * @see RectangularSystem
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
 * @brief Owns a constant right-handed orthonormal Cartesian orientation.
 *
 * The x direction anchors the frame. With two or three supplied directions,
 * Gram-Schmidt projection retains the transverse part of y and constructs z as
 * x cross y. A zero or numerically unusable transverse y is replaced by an
 * arbitrary orthogonal unit direction. The supplied z in the three-axis overload
 * does not determine the final triad because orthogonalize() recomputes it.
 * The one-axis overload completes y with Eigen's unitOrthogonal() operation.
 *
 * Input directions must be finite and x must be nonzero; the implementation has
 * no validation or fallback for an invalid x direction. Constructed axes and
 * cached transformations are private persistent definition data. For basis B,
 * local_to_global_ = B and global_to_local_ = B.transpose(). Point mappings are
 * pure rotations about the global origin, and get_axes() is position-independent.
 *
 * Euler helpers form Rx * Ry * Rz, which applies Rz first to a column vector,
 * then Ry and finally Rx. All angles and derivatives use radians. A transformed
 * copy rotates the triad, retains the name and ignores translation. Evaluation
 * and copying leave the source definition unchanged and own no solver state.
 */
class RectangularSystem : public CoordinateSystem {
public:
    // Build a constant triad; the final z direction is recomputed from x and y.
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis);
    RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis);
    RectangularSystem(const std::string& name, const Vec3& x_axis);

    // Pure point/component rotations about the global origin.
    Vec3 to_local(const Vec3& global_point) const override;
    Vec3 to_global(const Vec3& local_point) const override;
    // Constant global basis columns; the local evaluation point is ignored.
    Basis get_axes(const Vec3& local_point) const override;

    // Create an independent rotated orientation; translation is ignored.
    Ptr transformed(const Mat3& rotation, const Vec3& translation) const override;

    // Compose Rx(rot_x) * Ry(rot_y) * Rz(rot_z), with all angles in radians.
    static RectangularSystem euler(Precision rot_x, Precision rot_y, Precision rot_z);

    // Right-handed single-axis rotation matrices acting on column vectors.
    static Basis rotation_x(Precision angle);
    static Basis rotation_y(Precision angle);
    static Basis rotation_z(Precision angle);

    // Analytical derivatives of the single-axis matrices with respect to angle.
    static Basis rotation_x_derivative(Precision angle);
    static Basis rotation_y_derivative(Precision angle);
    static Basis rotation_z_derivative(Precision angle);

    // Partial derivatives of the complete Euler product with respect to each angle.
    static Basis derivative_rot_x(Precision rot_x, Precision rot_y, Precision rot_z);
    static Basis derivative_rot_y(Precision rot_x, Precision rot_y, Precision rot_z);
    static Basis derivative_rot_z(Precision rot_x, Precision rot_y, Precision rot_z);

private:
    // Prepare an orthonormal triad and cache its forward/inverse transformations.
    void normalize();
    void orthogonalize();
    void compute_transformations();

    // Cached component maps B and B.transpose(), with local axes in B columns.
    Basis local_to_global_{};
    Basis global_to_local_{};
    // Persistent right-handed unit directions expressed in global coordinates.
    Vec3 x_axis_{};
    Vec3 y_axis_{};
    Vec3 z_axis_{};
};
} // namespace cos
} // namespace fem
