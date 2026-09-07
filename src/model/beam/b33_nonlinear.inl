/**
 * @file b33_nonlinear.inl
 * @brief Implements finite-rotation nonlinear equilibrium for B33.
 *
 * The nonlinear B33 formulation keeps the existing two-node Euler-Bernoulli
 * element and its six nodal degrees of freedom. Translational coordinates are
 * evaluated in the current configuration and the three rotational coordinates
 * are interpreted as total global axis-angle vectors, identical to the finite-
 * rotation shell convention.
 *
 * A current corotated principal-section frame removes the complete rigid-body
 * motion of the element. The remaining axial extension, torsion and two bending
 * end rotations define one objective discrete strain energy whose infinitesimal
 * linearization is the existing B33 Euler-Bernoulli stiffness. Section reference
 * and shear-point offsets are rotated exactly instead of using the small-angle
 * rigid-offset matrices employed by the linear operator.
 *
 * Internal force and tangent are obtained as the first and second derivatives of
 * the same scalar strain energy. A compact second-order forward-differentiation
 * scalar is used locally for the twelve element degrees of freedom. The nodal
 * SO(3) derivatives themselves are not duplicated here: they are injected from
 * math::so3, so B33 and FRTShell use the same total-rotation convention and the
 * same exact first and second rotation derivatives.
 *
 * Linear stiffness, mass and classical prestress stiffness remain implemented
 * by the existing B33 paths in b33.h. This file is used only by nonlinear
 * equilibrium through B33::stiffness_tangent().
 *
 * @see B33
 * @see math::so3
 *
 * @author Finn Eggers
 * @date 07.09.2026
 */

#include "../../math/so3.h"

#include <array>
#include <cmath>

namespace fem::model {

namespace b33_nonlinear_detail {

// Number of independent element coordinates carried by the two-node B33.
static constexpr Index num_dofs = 12;

using Grad12 = StaticVector<num_dofs>;
using Hess12 = StaticMatrix<num_dofs, num_dofs>;

/**
 * @brief Scalar value with first and optional second derivatives for B33.
 *
 * The object carries derivatives with respect to the complete global element
 * coordinate vector
 *
 *     q = [u1x, u1y, u1z, r1x, r1y, r1z,
 *          u2x, u2y, u2z, r2x, r2y, r2z].
 *
 * Second derivatives are evaluated only when the caller requested a tangent.
 * Residual-only calls therefore retain the exact first derivative while avoiding
 * all Hessian products. This is deliberately local to B33: it is not intended as
 * a general automatic-differentiation layer for FEMaster.
 */
struct ADScalar {
    Precision value = Precision(0);
    Grad12    grad  = Grad12::Zero();
    Hess12    hess  = Hess12::Zero();
    bool      second_order = false;

    ADScalar() = default;

    explicit ADScalar(Precision value_in, bool second_order_in = false)
        : value(value_in), second_order(second_order_in) {}

    static ADScalar variable(Precision value, Index dof, bool second_order) {
        ADScalar result(value, second_order);
        result.grad(dof) = Precision(1);
        return result;
    }
};

// -----------------------------------------------------------------------------
// Scalar differentiation algebra
// -----------------------------------------------------------------------------

inline ADScalar operator+(const ADScalar& a, const ADScalar& b) {
    ADScalar result(a.value + b.value, a.second_order || b.second_order);
    result.grad = a.grad + b.grad;

    if (result.second_order) {
        result.hess = a.hess + b.hess;
    }
    return result;
}

inline ADScalar operator-(const ADScalar& a, const ADScalar& b) {
    ADScalar result(a.value - b.value, a.second_order || b.second_order);
    result.grad = a.grad - b.grad;

    if (result.second_order) {
        result.hess = a.hess - b.hess;
    }
    return result;
}

inline ADScalar operator-(const ADScalar& a) {
    ADScalar result(-a.value, a.second_order);
    result.grad = -a.grad;

    if (result.second_order) {
        result.hess = -a.hess;
    }
    return result;
}

inline ADScalar operator*(const ADScalar& a, const ADScalar& b) {
    ADScalar result(a.value * b.value, a.second_order || b.second_order);
    result.grad = b.value * a.grad + a.value * b.grad;

    if (result.second_order) {
        result.hess = b.value * a.hess
                    + a.value * b.hess
                    + a.grad * b.grad.transpose()
                    + b.grad * a.grad.transpose();
    }
    return result;
}

inline ADScalar operator*(Precision a, const ADScalar& b) {
    ADScalar result(a * b.value, b.second_order);
    result.grad = a * b.grad;

    if (result.second_order) {
        result.hess = a * b.hess;
    }
    return result;
}

inline ADScalar operator*(const ADScalar& a, Precision b) {
    return b * a;
}

inline ADScalar operator+(Precision a, const ADScalar& b) {
    ADScalar result = b;
    result.value += a;
    return result;
}

inline ADScalar operator+(const ADScalar& a, Precision b) {
    return b + a;
}

inline ADScalar operator-(const ADScalar& a, Precision b) {
    ADScalar result = a;
    result.value -= b;
    return result;
}

inline ADScalar operator-(Precision a, const ADScalar& b) {
    return ADScalar(a, b.second_order) - b;
}

/**
 * Applies one scalar function y = f(x) from its value and first two derivatives.
 *
 * For a scalar composition the multivariate chain rule is
 *
 *     grad(y) = f'(x) grad(x)
 *     Hess(y) = f'(x) Hess(x) + f''(x) grad(x) grad(x)^T.
 */
inline ADScalar compose(const ADScalar& x,
                        Precision       value,
                        Precision       first,
                        Precision       second) {
    ADScalar result(value, x.second_order);
    result.grad = first * x.grad;

    if (result.second_order) {
        result.hess = first * x.hess + second * x.grad * x.grad.transpose();
    }
    return result;
}

inline ADScalar inverse(const ADScalar& x) {
    logging::error(std::abs(x.value) > std::numeric_limits<Precision>::epsilon(),
        "B33: division by a vanishing nonlinear kinematic scalar");

    const Precision inv  = Precision(1) / x.value;
    const Precision inv2 = inv * inv;
    const Precision inv3 = inv2 * inv;
    return compose(x, inv, -inv2, Precision(2) * inv3);
}

inline ADScalar operator/(const ADScalar& a, const ADScalar& b) {
    return a * inverse(b);
}

inline ADScalar operator/(const ADScalar& a, Precision b) {
    logging::error(std::abs(b) > std::numeric_limits<Precision>::epsilon(),
        "B33: division by a vanishing constant");
    return (Precision(1) / b) * a;
}

inline ADScalar sqrt(const ADScalar& x) {
    logging::error(x.value > Precision(0),
        "B33: square root requires a positive nonlinear kinematic scalar");

    const Precision root   = std::sqrt(x.value);
    const Precision first  = Precision(0.5) / root;
    const Precision second = -Precision(0.25) / (x.value * root);
    return compose(x, root, first, second);
}

/**
 * Evaluates atan2(y,x) with complete first and second derivatives.
 *
 * The branch value comes from std::atan2 while the derivatives use the smooth
 * local differential away from x = y = 0. The relative-rotation logarithm calls
 * this only for a non-singular rotation state.
 */
inline ADScalar atan2(const ADScalar& y, const ADScalar& x) {
    const Precision radius2 = x.value * x.value + y.value * y.value;
    logging::error(radius2 > std::numeric_limits<Precision>::epsilon(),
        "B33: atan2 received a singular nonlinear rotation state");

    const Precision radius4 = radius2 * radius2;
    const Precision fx      = -y.value / radius2;
    const Precision fy      =  x.value / radius2;
    const Precision fxx     =  Precision(2) * x.value * y.value / radius4;
    const Precision fyy     = -Precision(2) * x.value * y.value / radius4;
    const Precision fxy     = (y.value * y.value - x.value * x.value) / radius4;

    ADScalar result(std::atan2(y.value, x.value), x.second_order || y.second_order);
    result.grad = fx * x.grad + fy * y.grad;

    if (result.second_order) {
        result.hess = fx * x.hess
                    + fy * y.hess
                    + fxx * x.grad * x.grad.transpose()
                    + fyy * y.grad * y.grad.transpose()
                    + fxy * (x.grad * y.grad.transpose() + y.grad * x.grad.transpose());
    }
    return result;
}

// -----------------------------------------------------------------------------
// Small three-dimensional vector and matrix operations for differentiated data
// -----------------------------------------------------------------------------

using ADVec3 = std::array<ADScalar, 3>;
using ADMat3 = std::array<std::array<ADScalar, 3>, 3>;

inline ADVec3 operator+(const ADVec3& a, const ADVec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline ADVec3 operator-(const ADVec3& a, const ADVec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline ADVec3 operator*(const ADScalar& a, const ADVec3& b) {
    return {a * b[0], a * b[1], a * b[2]};
}

inline ADVec3 operator*(Precision a, const ADVec3& b) {
    return {a * b[0], a * b[1], a * b[2]};
}

inline ADScalar dot(const ADVec3& a, const ADVec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline ADVec3 cross(const ADVec3& a, const ADVec3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

inline ADVec3 normalized(const ADVec3& value) {
    const ADScalar length = sqrt(dot(value, value));
    return {value[0] / length, value[1] / length, value[2] / length};
}

inline ADVec3 column(const ADMat3& matrix, Index col) {
    return {matrix[0][col], matrix[1][col], matrix[2][col]};
}

inline void set_column(ADMat3& matrix, Index col, const ADVec3& value) {
    matrix[0][col] = value[0];
    matrix[1][col] = value[1];
    matrix[2][col] = value[2];
}

/**
 * Multiplies a differentiated rotation by a constant global vector.
 */
inline ADVec3 multiply(const ADMat3& matrix, const Vec3& vector) {
    ADVec3 result;
    for (Index row = 0; row < 3; ++row) {
        result[row] = matrix[row][0] * vector(0)
                    + matrix[row][1] * vector(1)
                    + matrix[row][2] * vector(2);
    }
    return result;
}

/**
 * Multiplies a differentiated rotation by one constant reference basis.
 */
inline ADMat3 multiply(const ADMat3& matrix, const Mat3& basis) {
    ADMat3 result;
    for (Index row = 0; row < 3; ++row) {
        for (Index col = 0; col < 3; ++col) {
            result[row][col] = matrix[row][0] * basis(0, col)
                             + matrix[row][1] * basis(1, col)
                             + matrix[row][2] * basis(2, col);
        }
    }
    return result;
}

/**
 * Forms A^T * B for two differentiated orthonormal bases.
 */
inline ADMat3 transpose_multiply(const ADMat3& A, const ADMat3& B) {
    ADMat3 result;
    for (Index row = 0; row < 3; ++row) {
        for (Index col = 0; col < 3; ++col) {
            result[row][col] = A[0][row] * B[0][col]
                             + A[1][row] * B[1][col]
                             + A[2][row] * B[2][col];
        }
    }
    return result;
}

/**
 * Converts the existing analytical SO(3) derivatives into differentiated matrix
 * entries with respect to the complete B33 coordinate vector.
 *
 * The shell and beam therefore share exactly the same map
 *
 *     R(theta) = exp(skew(theta))
 *
 * and the same first and second derivatives. Only the placement into the twelve
 * beam coordinates is performed here.
 */
inline ADMat3 rotation_matrix(const Vec3& theta, Index first_dof, bool second_order) {
    Mat3                 rotation;
    std::array<Mat3, 3>  first;
    std::array<std::array<Mat3, 3>, 3> second;

    if (second_order) {
        math::so3::rotation_matrix_second_derivatives(theta, rotation, first, second);
    } else {
        math::so3::rotation_matrix_first_derivatives(theta, rotation, first);
    }

    ADMat3 result;
    for (Index row = 0; row < 3; ++row) {
        for (Index col = 0; col < 3; ++col) {
            ADScalar entry(rotation(row, col), second_order);

            for (Index a = 0; a < 3; ++a) {
                entry.grad(first_dof + a) = first[a](row, col);
            }

            if (second_order) {
                for (Index a = 0; a < 3; ++a) {
                    for (Index b = 0; b < 3; ++b) {
                        entry.hess(first_dof + a, first_dof + b) = second[a][b](row, col);
                    }
                }
            }
            result[row][col] = entry;
        }
    }
    return result;
}

/**
 * Evaluates the principal rotation vector of one relative rotation matrix.
 *
 * Let
 *
 *     v = 0.5 * axial(D - D^T) = sin(angle) * axis
 *     c = 0.5 * (trace(D) - 1) = cos(angle).
 *
 * Away from zero the rotation vector is
 *
 *     phi = atan2(|v|, c) / |v| * v.
 *
 * Around the identity the direct quotient is numerically singular although the
 * limit is smooth. The series
 *
 *     angle / sin(angle) = 1 + sin(angle)^2 / 6
 *                            + 3 sin(angle)^4 / 40 + ...
 *
 * is used there. A relative rotation arbitrarily close to 180 degrees has no
 * unique logarithm axis and is rejected explicitly. Large rigid-body rotations
 * do not approach this singularity because they are removed by the corotated
 * frame before the relative rotations are formed.
 */
inline ADVec3 rotation_log(const ADMat3& rotation) {
    ADVec3 axial = {
        Precision(0.5) * (rotation[2][1] - rotation[1][2]),
        Precision(0.5) * (rotation[0][2] - rotation[2][0]),
        Precision(0.5) * (rotation[1][0] - rotation[0][1])
    };

    const ADScalar sin_squared = dot(axial, axial);
    const ADScalar cosine = Precision(0.5)
        * (rotation[0][0] + rotation[1][1] + rotation[2][2] - Precision(1));

    ADScalar scale;

    // Around the identity use the regular series directly in sin(angle)^2.
    if (sin_squared.value < Precision(1e-10) && cosine.value > Precision(0)) {
        scale = Precision(1)
              + sin_squared / Precision(6)
              + Precision(3) * sin_squared * sin_squared / Precision(40);
    } else {
        logging::error(sin_squared.value > Precision(1e-14),
            "B33: relative section rotation is too close to 180 degrees");

        const ADScalar sine  = sqrt(sin_squared);
        const ADScalar angle = atan2(sine, cosine);
        scale = angle / sine;
    }

    return scale * axial;
}

} // namespace b33_nonlinear_detail

/**
 * Evaluates the objective finite-rotation B33 internal force and tangent.
 *
 * The nonlinear element is expressed through one scalar strain energy. The
 * construction proceeds in the following order:
 *
 * 1. Build exact current nodal section rotations from the total axis-angle DOFs.
 * 2. Rotate the reference-to-shear-point offset with each nodal section.
 * 3. Construct the current shear-point chord and one objective corotated frame.
 * 4. Remove that rigid-body frame from both nodal section orientations.
 * 5. Measure axial extension and the two finite relative section rotations.
 * 6. Insert those deformation measures into the classical B33 basic energy.
 * 7. Differentiate the same energy once for internal force and twice for tangent.
 *
 * The corotated frame uses the current chord as its local x axis. Its local y
 * direction is the projected average of both current principal y axes. Thus a
 * common finite rigid rotation, including a common twist about the beam axis,
 * rotates the frame and both nodal sections identically and produces exactly
 * zero deformation energy.
 *
 * The basic bending energy is written in terms of the two end rotations relative
 * to the current chord. In the infinitesimal limit these become
 *
 *     phi1_z = theta1_z - (v2 - v1) / L
 *     phi2_z = theta2_z - (v2 - v1) / L
 *
 * and
 *
 *     phi1_y = theta1_y + (w2 - w1) / L
 *     phi2_y = theta2_y + (w2 - w1) / L.
 *
 * Applying the rotational B33 submatrix [4,2;2,4] * EI/L to these four basic
 * rotations reproduces the complete cubic Euler-Bernoulli translational and
 * rotational stiffness after linearization. The same construction is used for
 * both principal bending planes. Torsion depends only on phi2_x - phi1_x, and
 * axial response depends on the exact current shear-point chord extension.
 *
 * Section offsets are treated geometrically. If r_ref_sp denotes the constant
 * reference vector from the nodal reference point to the shear point, the
 * current mechanical point is
 *
 *     x_sp_i = X_ref_i + u_i + R_i * r_ref_sp.
 *
 * Its infinitesimal variation is the same rigid-offset relation used by the
 * existing linear B33 matrix, while finite rotations remain exact.
 *
 * @param buffer Optional caller-owned 12 x 12 tangent storage. A null pointer
 *               requests residual only.
 * @param nodal_forces Global nodal internal-force accumulator.
 * @param displacement Trial nodal displacement and total-rotation field.
 * @return Mapped consistent tangent, or an empty map for residual-only calls.
 */
inline MapMatrix B33::stiffness_tangent(Precision*   buffer,
                                        NodeData&    nodal_forces,
                                        const Field& displacement) {
    using namespace b33_nonlinear_detail;

    // -------------------------------------------------------------------------
    // Validate the nonlinear element state and gather the twelve trial DOFs
    // -------------------------------------------------------------------------

    logging::error(displacement.domain == FieldDomain::NODE,
        "B33: nonlinear displacement field must use NODE domain");
    logging::error(displacement.components >= 6,
        "B33: nonlinear displacement field requires six nodal components");
    logging::error(nodal_forces.components >= 6,
        "B33: nonlinear internal force requires six nodal components");

    const bool with_tangent = buffer != nullptr;

    StaticVector<12> q;
    for (Index node = 0; node < 2; ++node) {
        q.template segment<6>(6 * node) =
            displacement.row_vec6(static_cast<Index>(node_ids[node]));
    }

    // -------------------------------------------------------------------------
    // Build the reference principal-section geometry used by the existing B33
    // -------------------------------------------------------------------------

    const Vec3 X1 = this->node_position(0);
    const Vec3 X2 = this->node_position(1);
    const Precision L = (X2 - X1).norm();

    logging::error(L > std::numeric_limits<Precision>::epsilon(),
        "B33: nonlinear beam requires a positive reference length for element ", elem_id);

    Profile* profile = get_profile();

    const Precision E  = get_elasticity()->youngs;
    const Precision G  = get_elasticity()->shear;
    const Precision A  = profile->area_;
    const Precision It = profile->torsion_inertia_;

    Precision Iy  = profile->inertia_y_;
    Precision Iz  = profile->inertia_z_;
    Precision Iyz = profile->product_inertia_yz_;

    // The linear B33 diagonalizes the y-z inertia tensor before constructing the
    // uncoupled bending blocks. Use the identical principal angle and principal
    // inertias here so the nonlinear energy has exactly the same small-strain
    // material operator.
    const Precision principal_phi = principal_angle();
    const Precision inertia_scale = std::max<Precision>(
        Precision(1), std::abs(Iy) + std::abs(Iz));

    if (std::abs(Iyz) > inertia_scale * Precision(1e-14)) {
        const Precision c  = std::cos(principal_phi);
        const Precision s  = std::sin(principal_phi);
        const Precision c2 = c * c;
        const Precision s2 = s * s;
        const Precision sc = s * c;

        const Precision Iy_principal = Iy * c2 + Iz * s2 - Precision(2) * Iyz * sc;
        const Precision Iz_principal = Iy * s2 + Iz * c2 + Precision(2) * Iyz * sc;
        Iy = Iy_principal;
        Iz = Iz_principal;
    }

    // principal_rotation_matrix() maps global vectors into local principal
    // coordinates. Its transpose therefore stores the three reference principal
    // section axes as global columns [e1, e2, e3].
    const Mat3 reference_basis = principal_rotation_matrix().transpose();

    // The linear element first maps REF -> SMP and then SMP -> SP. Because both
    // are rigid translations in the same section frame their composition is one
    // vector REF -> SP. Rotate the stored offsets into the same principal frame
    // before forming that vector.
    Precision ey   = profile->offset_y_;
    Precision ez   = profile->offset_z_;
    Precision refy = profile->reference_y_;
    Precision refz = profile->reference_z_;

    BeamElement<2>::rotate_yz_to_principal(principal_phi, ey, ez);
    BeamElement<2>::rotate_yz_to_principal(principal_phi, refy, refz);

    Vec3 offset_local;
    offset_local << Precision(0), ey - refy, ez - refz;
    const Vec3 offset_global = reference_basis * offset_local;

    // -------------------------------------------------------------------------
    // Construct differentiated current nodal positions and section orientations
    // -------------------------------------------------------------------------

    ADVec3 x_ref_1;
    ADVec3 x_ref_2;
    for (Index component = 0; component < 3; ++component) {
        x_ref_1[component] = ADScalar::variable(
            X1(component) + q(component), component, with_tangent);
        x_ref_2[component] = ADScalar::variable(
            X2(component) + q(6 + component), 6 + component, with_tangent);
    }

    const Vec3 theta_1 = q.template segment<3>(3);
    const Vec3 theta_2 = q.template segment<3>(9);

    // Reuse the exact SO(3) first and second derivatives already employed by the
    // finite-rotation shells. Rotational coordinates are therefore interpreted
    // identically across shell and beam nonlinear formulations.
    const ADMat3 R1 = b33_nonlinear_detail::rotation_matrix(theta_1, 3, with_tangent);
    const ADMat3 R2 = b33_nonlinear_detail::rotation_matrix(theta_2, 9, with_tangent);

    const ADMat3 basis_1 = multiply(R1, reference_basis);
    const ADMat3 basis_2 = multiply(R2, reference_basis);

    // Rotate the complete reference-to-shear-point vector with the nodal section.
    // This is the finite-rotation counterpart of the linear rigid_offset_N()
    // transformations used by stiffness_impl().
    const ADVec3 x_sp_1 = x_ref_1 + multiply(R1, offset_global);
    const ADVec3 x_sp_2 = x_ref_2 + multiply(R2, offset_global);

    // -------------------------------------------------------------------------
    // Build the objective current corotated principal-section frame
    // -------------------------------------------------------------------------

    const ADVec3 chord = x_sp_2 - x_sp_1;
    const ADScalar current_length = sqrt(dot(chord, chord));

    logging::error(current_length.value > std::numeric_limits<Precision>::epsilon(),
        "B33: nonlinear current beam length vanished for element ", elem_id);

    const ADVec3 e1 = {
        chord[0] / current_length,
        chord[1] / current_length,
        chord[2] / current_length
    };

    // The chord fixes only two rotational components. A third condition is needed
    // to remove rigid twist about the chord. Use the average current principal y
    // direction and project it into the plane normal to e1. Under any common
    // rigid-body rotation both the chord and the nodal y axes rotate by the same
    // Q, so the resulting corotated frame also transforms exactly by Q.
    ADVec3 y_seed = column(basis_1, 1) + column(basis_2, 1);
    ADVec3 y_projected = y_seed - dot(e1, y_seed) * e1;
    ADScalar y_projected_norm2 = dot(y_projected, y_projected);

    // Opposite nodal y axes make their average undefined. This is the discrete
    // mean-frame singularity associated with approximately 180 degrees of
    // relative twist. Use node 1 as a deterministic local fallback so the frame
    // remains defined up to the separate relative-logarithm singularity check.
    if (y_projected_norm2.value < Precision(1e-12)) {
        y_seed = column(basis_1, 1);
        y_projected = y_seed - dot(e1, y_seed) * e1;
        y_projected_norm2 = dot(y_projected, y_projected);
    }

    logging::error(y_projected_norm2.value > Precision(1e-14),
        "B33: unable to construct nonlinear corotated section frame for element ", elem_id);

    ADVec3 e2 = normalized(y_projected);
    ADVec3 e3 = normalized(cross(e1, e2));

    // Recompute e2 from e3 x e1. This removes the last round-off component along
    // e1 and gives one right-handed orthonormal frame [e1,e2,e3].
    e2 = cross(e3, e1);

    ADMat3 corotated_basis;
    set_column(corotated_basis, 0, e1);
    set_column(corotated_basis, 1, e2);
    set_column(corotated_basis, 2, e3);

    // -------------------------------------------------------------------------
    // Remove rigid-body rotation and obtain finite basic deformation measures
    // -------------------------------------------------------------------------

    const ADMat3 relative_1 = transpose_multiply(corotated_basis, basis_1);
    const ADMat3 relative_2 = transpose_multiply(corotated_basis, basis_2);

    const ADVec3 phi_1 = rotation_log(relative_1);
    const ADVec3 phi_2 = rotation_log(relative_2);

    // Exact axial extension of the mechanical shear-point line. The reference
    // offset is identical at both nodes, so its reference length is the nodal B33
    // reference length L. Rotated offsets nevertheless contribute exactly in the
    // current configuration.
    const ADScalar extension = current_length - L;

    // -------------------------------------------------------------------------
    // Assemble the objective discrete Euler-Bernoulli strain energy
    // -------------------------------------------------------------------------

    // Axial energy:
    //     Pi_axial = 1/2 * EA/L * extension^2.
    ADScalar energy = Precision(0.5) * (E * A / L) * extension * extension;

    // Torsional energy depends only on the relative twist between both sections.
    // A common twist is absorbed by the corotated frame and therefore contributes
    // exactly zero energy.
    const ADScalar twist = phi_2[0] - phi_1[0];
    energy = energy + Precision(0.5) * (G * It / L) * twist * twist;

    // Euler-Bernoulli bending in the principal y plane. The two end rotations
    // relative to the chord are acted on by the rotational basic stiffness
    //
    //     EI/L * [4 2; 2 4].
    //
    // Writing the quadratic form directly avoids creating a temporary matrix.
    energy = energy + Precision(2) * (E * Iy / L)
        * (phi_1[1] * phi_1[1]
         + phi_1[1] * phi_2[1]
         + phi_2[1] * phi_2[1]);

    // Identical bending energy in the principal z plane.
    energy = energy + Precision(2) * (E * Iz / L)
        * (phi_1[2] * phi_1[2]
         + phi_1[2] * phi_2[2]
         + phi_2[2] * phi_2[2]);

    // -------------------------------------------------------------------------
    // Scatter the energy gradient as internal force and expose the Hessian
    // -------------------------------------------------------------------------

    // Because force and tangent are derivatives of exactly the same scalar
    // potential, the nonlinear residual is conservative and the tangent is
    // analytically consistent with it. No independently derived geometric
    // stiffness is added in this path; all stress-dependent terms are already
    // contained in the Hessian of the finite-rotation energy.
    const Grad12 internal_force = energy.grad;

    for (Index node = 0; node < 2; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);
        for (Index dof = 0; dof < 6; ++dof) {
            nodal_forces(node_id, dof) += internal_force(6 * node + dof);
        }
    }

    if (!with_tangent) {
        return MapMatrix(nullptr, 0, 0);
    }

    // The Hessian is symmetric by construction. Average the two triangles only
    // to remove floating-point asymmetry accumulated by the local matrix algebra.
    const Hess12 tangent = Precision(0.5) * (energy.hess + energy.hess.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = tangent;
    return mapped;
}

} // namespace fem::model
