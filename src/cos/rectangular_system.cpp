/**
 * @file rectangular_system.cpp
 * @brief Implements Cartesian basis construction and analytical rotations.
 *
 * The coordinate-system implementation completes right-handed orthonormal axes
 * with Gram-Schmidt projection and caches the basis matrix and its transpose.
 * Point/component mappings use these constant matrices without a spatial offset.
 * Rigid placement creates a rotated independent definition.
 *
 * Single-axis rotations and their analytical derivatives support the fixed Euler
 * product Rx * Ry * Rz. Consumers such as elements own the constitutive operations
 * and sensitivities that use these matrices; this file changes no solver state.
 *
 * @see RectangularSystem
 * @see CoordinateSystem
 *
 * @author Finn Eggers
 * @date 06.03.2025
 */

#include "rectangular_system.h"

#include <Eigen/Geometry>

namespace fem {
namespace cos {

/**
 * @brief Builds a Cartesian orientation from supplied global directions.
 *
 * Orthogonalize x/y, recompute z = x cross y, normalize the triad and cache B and
 * B.transpose(). The supplied z is initially stored but replaced by the cross
 * product; it does not select the final handedness or orientation. The x input
 * must be finite and nonzero. Degenerate transverse y uses the fallback described
 * by orthogonalize(); no validation of invalid x is performed.
 *
 * @param name Immutable definition identifier.
 * @param x_axis Global direction anchoring the local x axis.
 * @param y_axis Global direction whose transverse component anchors local y.
 * @param z_axis Initially supplied global z direction, recomputed during construction.
 */
RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis, const Vec3& z_axis)
    : CoordinateSystem(name)
    , x_axis_(x_axis)
    , y_axis_(y_axis)
    , z_axis_(z_axis) {
    // Complete a right-handed unit triad before caching the component maps.
    orthogonalize();
    normalize();
    compute_transformations();
}

/**
 * @brief Completes a Cartesian orientation from global x/y directions.
 *
 * Initialize z with x cross y, then orthogonalize and normalize the triad and
 * cache its forward/inverse component maps. The final z is recomputed from the
 * prepared x/y axes. Finite nonzero x is required; zero or unusable transverse y
 * is completed by an arbitrary orthogonal direction in orthogonalize().
 *
 * @param name Immutable definition identifier.
 * @param x_axis Global direction anchoring the local x axis.
 * @param y_axis Global direction defining the desired transverse y orientation.
 */
RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis, const Vec3& y_axis)
    : CoordinateSystem(name)
    , x_axis_(x_axis)
    , y_axis_(y_axis)
    , z_axis_(x_axis.cross(y_axis)) {
    // Complete a right-handed unit triad before caching the component maps.
    orthogonalize();
    normalize();
    compute_transformations();
}

/**
 * @brief Completes a Cartesian orientation from one global x direction.
 *
 * Normalize x, obtain an arbitrary orthogonal unit y from Eigen and form the
 * right-handed unit z by a cross product. Cache B and B.transpose() for subsequent
 * read-only evaluation. Input x must be finite and nonzero; no invalid-direction
 * check is performed and the transverse orientation is not user-prescribed.
 *
 * @param name Immutable definition identifier.
 * @param x_axis Finite nonzero global direction of the local x axis.
 */
RectangularSystem::RectangularSystem(const std::string& name, const Vec3& x_axis)
    : CoordinateSystem(name), x_axis_(x_axis) {
    // Complete the transverse plane from a valid x direction without prescribing roll.
    x_axis_.normalize();
    y_axis_ = x_axis_.unitOrthogonal();
    z_axis_ = x_axis_.cross(y_axis_).normalized();
    compute_transformations();
}

/**
 * @brief Resolves a global Cartesian point or vector into local components.
 *
 * For the constant orthonormal basis B, apply B.transpose(). There is no origin
 * subtraction, so positions rotate about the global origin just like vectors.
 * The definition and cached matrices remain unchanged.
 *
 * @param global_point Global Cartesian components.
 * @return Local Cartesian components in the stored x/y/z basis.
 */
Vec3 RectangularSystem::to_local(const Vec3& global_point) const {
    return global_to_local_ * global_point;
}

/**
 * @brief Maps local Cartesian components into the global frame.
 *
 * Apply the constant basis B with local unit directions stored as columns. There
 * is no translation. This operation inverts to_local() for a valid orthonormal
 * frame and leaves the definition unchanged.
 *
 * @param local_point Local Cartesian point or vector components.
 * @return Global Cartesian components.
 */
Vec3 RectangularSystem::to_global(const Vec3& local_point) const {
    return local_to_global_ * local_point;
}

/**
 * @brief Returns the position-independent Cartesian basis.
 *
 * @param local_point Unused; a rectangular orientation has constant axes.
 * @return Cached matrix B whose columns are local x/y/z unit directions in global
 *          coordinates. The definition remains unchanged.
 */
Basis RectangularSystem::get_axes(const Vec3& local_point) const {
    (void)local_point;
    return local_to_global_;
}

/**
 * @brief Creates an independently owned orientation under rigid placement.
 *
 * Rotate all stored unit directions and reconstruct the Cartesian definition.
 * Translation is ignored because the system has no spatial origin. The supplied
 * rotation must be proper orthonormal and finite; no validation is performed here.
 * Reconstruction orthogonalizes and normalizes the rotated triad.
 *
 * @param rotation Proper orthonormal placement rotation in global coordinates.
 * @param translation Unused instance translation.
 * @return Independent rotated definition retaining the name; source state is unchanged.
 */
CoordinateSystem::Ptr RectangularSystem::transformed(const Mat3& rotation, const Vec3& translation) const {
    // An orientation has no translated origin; rotate only its global unit axes.
    (void) translation;
    return std::make_shared<RectangularSystem>(
        name,
        rotation * x_axis_,
        rotation * y_axis_,
        rotation * z_axis_
    );
}

/**
 * @brief Constructs an unnamed orientation from the fixed Euler product.
 *
 * Form R = Rx(rot_x) * Ry(rot_y) * Rz(rot_z) and use its columns as global axes.
 * For column vectors, the rightmost z rotation acts first, then y and finally x.
 * The constructor prepares the orthonormal basis and caches both component maps.
 * No source object is modified; angles must be finite.
 *
 * @param rot_x Right-handed x rotation angle in radians.
 * @param rot_y Right-handed y rotation angle in radians.
 * @param rot_z Right-handed z rotation angle in radians.
 * @return Cartesian definition with an empty name and basis R.
 */
RectangularSystem RectangularSystem::euler(Precision rot_x, Precision rot_y, Precision rot_z) {
    // Compose the fixed rotation order before extracting the local basis columns.
    Basis rot_mat = rotation_x(rot_x) * rotation_y(rot_y) * rotation_z(rot_z);
    return RectangularSystem("", rot_mat.col(0), rot_mat.col(1), rot_mat.col(2));
}

/**
 * @brief Returns the right-handed rotation matrix about the x axis.
 *
 * The standard sine/cosine matrix acts on Cartesian column vectors and leaves the
 * x component unchanged. Its columns give the rotated Cartesian unit directions.
 * The angle must be finite; the helper performs no validation or state update.
 *
 * @param angle Rotation angle in radians.
 * @return Proper orthonormal matrix Rx(angle).
 */
Basis RectangularSystem::rotation_x(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << 1, 0, 0,
           0, std::cos(angle), -std::sin(angle),
           0, std::sin(angle), std::cos(angle);
    return rot;
}

/**
 * @brief Returns the right-handed rotation matrix about the y axis.
 *
 * The standard sine/cosine matrix acts on Cartesian column vectors and leaves the
 * y component unchanged. Its columns give the rotated Cartesian unit directions.
 * The angle must be finite; the helper performs no validation or state update.
 *
 * @param angle Rotation angle in radians.
 * @return Proper orthonormal matrix Ry(angle).
 */
Basis RectangularSystem::rotation_y(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << std::cos(angle), 0, std::sin(angle),
           0, 1, 0,
           -std::sin(angle), 0, std::cos(angle);
    return rot;
}

/**
 * @brief Returns the right-handed rotation matrix about the z axis.
 *
 * The standard sine/cosine matrix acts on Cartesian column vectors and leaves the
 * z component unchanged. Its columns give the rotated Cartesian unit directions.
 * The angle must be finite; the helper performs no validation or state update.
 *
 * @param angle Rotation angle in radians.
 * @return Proper orthonormal matrix Rz(angle).
 */
Basis RectangularSystem::rotation_z(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << std::cos(angle), -std::sin(angle), 0,
           std::sin(angle), std::cos(angle), 0,
           0, 0, 1;
    return rot;
}

/**
 * @brief Differentiates the single-axis x rotation with respect to its angle.
 *
 * Differentiate the sine/cosine entries of Rx(angle); constant entries have zero
 * derivative. The result maps column vectors to their angular rate per radian and
 * is not itself a rotation matrix. No definition state is modified.
 *
 * @param angle Finite rotation angle in radians.
 * @return Analytical matrix derivative dRx/dangle.
 */
Basis RectangularSystem::rotation_x_derivative(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << 0, 0, 0,
           0, -std::sin(angle), -std::cos(angle),
           0, std::cos(angle), -std::sin(angle);
    return rot;
}

/**
 * @brief Differentiates the single-axis y rotation with respect to its angle.
 *
 * Differentiate the sine/cosine entries of Ry(angle); constant entries have zero
 * derivative. The result maps column vectors to their angular rate per radian and
 * is not itself a rotation matrix. No definition state is modified.
 *
 * @param angle Finite rotation angle in radians.
 * @return Analytical matrix derivative dRy/dangle.
 */
Basis RectangularSystem::rotation_y_derivative(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << -std::sin(angle), 0, std::cos(angle),
           0, 0, 0,
           -std::cos(angle), 0, -std::sin(angle);
    return rot;
}

/**
 * @brief Differentiates the single-axis z rotation with respect to its angle.
 *
 * Differentiate the sine/cosine entries of Rz(angle); constant entries have zero
 * derivative. The result maps column vectors to their angular rate per radian and
 * is not itself a rotation matrix. No definition state is modified.
 *
 * @param angle Finite rotation angle in radians.
 * @return Analytical matrix derivative dRz/dangle.
 */
Basis RectangularSystem::rotation_z_derivative(Precision angle) {
    // Assemble the analytical Cartesian matrix in row/column component order.
    Basis rot;
    rot << -std::sin(angle), -std::cos(angle), 0,
           std::cos(angle), -std::sin(angle), 0,
           0, 0, 0;
    return rot;
}

/**
 * @brief Computes the Euler-product partial derivative with respect to rot_x.
 *
 * For R = Rx(rot_x) * Ry(rot_y) * Rz(rot_z), differentiate only the x factor:
 * dR/drot_x = Rx' * Ry * Rz. Preserve the order because rotations
 * and their derivatives do not generally commute. All angles must be finite;
 * this analytical sensitivity changes no coordinate-system state.
 *
 * @param rot_x Right-handed x rotation angle in radians.
 * @param rot_y Right-handed y rotation angle in radians.
 * @param rot_z Right-handed z rotation angle in radians.
 * @return Partial derivative of the complete basis matrix per radian.
 */
Basis RectangularSystem::derivative_rot_x(Precision rot_x, Precision rot_y, Precision rot_z) {
    // Differentiate the selected factor while retaining the Euler composition order.
    return rotation_x_derivative(rot_x) * rotation_y(rot_y) * rotation_z(rot_z);
}

/**
 * @brief Computes the Euler-product partial derivative with respect to rot_y.
 *
 * For R = Rx(rot_x) * Ry(rot_y) * Rz(rot_z), differentiate only the y factor:
 * dR/drot_y = Rx * Ry' * Rz. Preserve the order because rotations
 * and their derivatives do not generally commute. All angles must be finite;
 * this analytical sensitivity changes no coordinate-system state.
 *
 * @param rot_x Right-handed x rotation angle in radians.
 * @param rot_y Right-handed y rotation angle in radians.
 * @param rot_z Right-handed z rotation angle in radians.
 * @return Partial derivative of the complete basis matrix per radian.
 */
Basis RectangularSystem::derivative_rot_y(Precision rot_x, Precision rot_y, Precision rot_z) {
    // Differentiate the selected factor while retaining the Euler composition order.
    return rotation_x(rot_x) * rotation_y_derivative(rot_y) * rotation_z(rot_z);
}

/**
 * @brief Computes the Euler-product partial derivative with respect to rot_z.
 *
 * For R = Rx(rot_x) * Ry(rot_y) * Rz(rot_z), differentiate only the z factor:
 * dR/drot_z = Rx * Ry * Rz'. Preserve the order because rotations
 * and their derivatives do not generally commute. All angles must be finite;
 * this analytical sensitivity changes no coordinate-system state.
 *
 * @param rot_x Right-handed x rotation angle in radians.
 * @param rot_y Right-handed y rotation angle in radians.
 * @param rot_z Right-handed z rotation angle in radians.
 * @return Partial derivative of the complete basis matrix per radian.
 */
Basis RectangularSystem::derivative_rot_z(Precision rot_x, Precision rot_y, Precision rot_z) {
    // Differentiate the selected factor while retaining the Euler composition order.
    return rotation_x(rot_x) * rotation_y(rot_y) * rotation_z_derivative(rot_z);
}

/**
 * @brief Normalizes the stored Cartesian axes in place.
 *
 * Each direction is divided by its Euclidean norm. Orthogonality and handedness
 * must already be established by orthogonalize(). Zero and nonfinite directions
 * are not checked. Cached matrices are not updated here; construction subsequently
 * calls compute_transformations() to make them consistent with the prepared axes.
 */
void RectangularSystem::normalize() {
    // Enforce unit lengths after the triad has been orthogonalized.
    x_axis_.normalize();
    y_axis_.normalize();
    z_axis_.normalize();
}

/**
 * @brief Prepares a right-handed triad with x as the fixed reference direction.
 *
 * Normalize x and Gram-Schmidt project y with y_perp = y - (x dot y) * x.
 * If y is zero, use an arbitrary unit orthogonal direction. Otherwise normalize
 * y_perp and replace a nonfinite or negligibly small result by the same fallback.
 * Finally recompute z = normalize(x cross y), replacing any supplied z direction.
 *
 * A finite nonzero x is a caller precondition and has no fallback. This method
 * updates the stored axes only; cached component maps are built separately.
 */
void RectangularSystem::orthogonalize() {
    // Anchor the frame to unit x and choose an orthogonal y if no direction was supplied.
    x_axis_.normalize();
    if (y_axis_.norm() == 0) y_axis_ = x_axis_.unitOrthogonal();
    else {
        // Remove the component along x; a collinear y can become nonfinite on normalization.
        y_axis_ = (y_axis_ - (x_axis_.dot(y_axis_) * x_axis_)).normalized();
        if (!std::isfinite(y_axis_.squaredNorm()) || y_axis_.squaredNorm() < 1e-30)
            y_axis_ = x_axis_.unitOrthogonal();
    }
    // Complete a right-handed unit z, independent of its initial supplied direction.
    z_axis_ = x_axis_.cross(y_axis_).normalized();
}

/**
 * @brief Caches both Cartesian component maps from the prepared unit axes.
 *
 * Store the global unit directions as columns of B = local_to_global_. For an
 * orthonormal triad B inverse equals B.transpose(), so global_to_local_ is obtained
 * without a matrix inversion. The axes must already be normalized and orthogonal;
 * this method only updates the two persistent transformation matrices.
 */
void RectangularSystem::compute_transformations() {
    // Build B from global axis columns and use orthonormality for its inverse.
    local_to_global_.col(0) = x_axis_;
    local_to_global_.col(1) = y_axis_;
    local_to_global_.col(2) = z_axis_;
    global_to_local_ = local_to_global_.transpose();
}
} // namespace cos
} // namespace fem
