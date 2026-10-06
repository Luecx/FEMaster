/**
 * @file b31.cpp
 * @brief Implements the two-node geometrically exact B31 beam.
 *
 * The formulation uses a reference-line description with exact finite rotations.
 * Translation/shear strains are evaluated at the section shear point and the
 * relative nodal rotation is interpolated on SO(3) along its shortest geodesic.
 * This keeps a superposed rigid-body motion exactly strain free.
 */

#include "b31.h"

#include "../../core/logging.h"
#include "../../math/so3.h"
#include "../../math/vec_util.h"
#include "../../section/profile.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <Eigen/Geometry>

namespace fem::model {

namespace {

Vec3 vee_skew(const Mat3& matrix) {
    return Vec3(
        Precision(0.5) * (matrix(2, 1) - matrix(1, 2)),
        Precision(0.5) * (matrix(0, 2) - matrix(2, 0)),
        Precision(0.5) * (matrix(1, 0) - matrix(0, 1))
    );
}

Vec3 rotation_vector_from_matrix(const Mat3& rotation) {
    Eigen::AngleAxis<Precision> angle_axis(rotation);
    if (std::abs(angle_axis.angle()) <= std::numeric_limits<Precision>::epsilon()) {
        return Vec3::Zero();
    }
    return angle_axis.angle() * angle_axis.axis();
}

Mat3 right_jacobian_inverse(const Vec3& rotation_vector) {
    const Precision angle_squared = rotation_vector.squaredNorm();
    const Mat3 K = math::skew(rotation_vector);

    if (angle_squared < Precision(1e-8)) {
        return Mat3::Identity()
             + Precision(0.5) * K
             + (Precision(1) / Precision(12)
                + angle_squared / Precision(720)) * K * K;
    }

    const Precision angle = std::sqrt(angle_squared);
    const Precision coefficient =
        Precision(1) / angle_squared
        - (Precision(1) + std::cos(angle))
          / (Precision(2) * angle * std::sin(angle));

    return Mat3::Identity()
         + Precision(0.5) * K
         + coefficient * K * K;
}


/**
 * Directional derivative of the inverse SO(3) right Jacobian.
 *
 * For phi = Log(R) and the body increment eta = vee(R^T dR),
 *
 *     d phi = J_r(phi)^-1 eta.
 *
 * The derivative below is required by the exact second derivative of the
 * relative-rotation logarithm used by the B31 curvature field.
 */
Mat3 right_jacobian_inverse_directional(
    const Vec3& rotation_vector,
    const Vec3& direction
) {
    const Precision angle_squared = rotation_vector.squaredNorm();
    const Precision projection    = rotation_vector.dot(direction);

    const Mat3 K  = math::skew(rotation_vector);
    const Mat3 dK = math::skew(direction);

    Precision coefficient            = Precision(0);
    Precision coefficient_directional = Precision(0);

    if (angle_squared < Precision(1e-8)) {
        const Precision s2 = angle_squared * angle_squared;

        coefficient =
              Precision(1) / Precision(12)
            + angle_squared / Precision(720)
            + s2 / Precision(30240);

        coefficient_directional =
            (Precision(1) / Precision(360)
             + angle_squared / Precision(7560)
             + s2 / Precision(201600)) * projection;
    } else {
        const Precision angle      = std::sqrt(angle_squared);
        const Precision half       = Precision(0.5) * angle;
        const Precision sin_half   = std::sin(half);
        const Precision cos_half   = std::cos(half);
        const Precision cot_half   = cos_half / sin_half;
        const Precision csc2_half  = Precision(1) / (sin_half * sin_half);

        coefficient =
              Precision(1) / angle_squared
            - cot_half / (Precision(2) * angle);

        const Precision dcoefficient_dangle =
              -Precision(2) / (angle * angle_squared)
              + csc2_half / (Precision(4) * angle)
              + cot_half / (Precision(2) * angle_squared);

        coefficient_directional =
            dcoefficient_dangle * projection / angle;
    }

    return Precision(0.5) * dK
         + coefficient_directional * K * K
         + coefficient * (dK * K + K * dK);
}

} // namespace

Vec3 B31::node_position_reference(Index local_node) const {
    logging::error(local_node < N,
        "B31: local node index out of range in element ", this->elem_id);
    logging::error(this->_model_data != nullptr,
        "B31: no model data assigned to element ", this->elem_id);
    logging::error(this->_model_data->positions_reference != nullptr,
        "B31: reference positions field is not initialized");

    return this->_model_data->positions_reference->row_vec3(
        static_cast<Index>(node_ids[local_node]));
}

Precision B31::reference_length() const {
    return (node_position_reference(1) - node_position_reference(0)).norm();
}

Mat3 B31::reference_frame() {
    const Vec3 axis = node_position_reference(1) - node_position_reference(0);
    const Precision length = axis.norm();
    logging::error(length > Precision(0),
        "B31: zero reference length in element ", this->elem_id);

    const Vec3 e1 = axis / length;
    const Vec3 n1 = orientation_direction();
    const Vec3 cross = e1.cross(n1);
    logging::error(cross.norm() > Precision(1e-12),
        "B31: beam orientation is parallel to the reference axis in element ",
        this->elem_id);

    const Vec3 e3 = cross.normalized();
    const Vec3 e2 = e3.cross(e1).normalized();

    Mat3 frame;
    frame.col(0) = e1;
    frame.col(1) = e2;
    frame.col(2) = e3;
    return frame;
}

B31::Vec12 B31::element_displacement(const Field* displacement) const {
    Vec12 q = Vec12::Zero();
    if (!displacement) return q;

    logging::error(displacement->domain == FieldDomain::NODE,
        "B31: displacement field must use NODE domain");
    logging::error(displacement->components >= dofs_per_node,
        "B31: displacement field requires six components");

    for (Index node = 0; node < N; ++node) {
        q.template segment<6>(node * dofs_per_node) =
            displacement->row_vec6(static_cast<Index>(node_ids[node]));
    }
    return q;
}

B31::Kinematics B31::kinematics(const Vec12& q, bool with_B) {
    const Precision L = reference_length();
    const Mat3 Qref = reference_frame();

    const Vec3 e1_local = Vec3::UnitX();

    const Profile* profile = get_profile();
    const Vec3 offset_ref_to_sp_local(
        Precision(0),
        profile->offset_y_ - profile->reference_y_,
        profile->offset_z_ - profile->reference_z_);
    const Vec3 offset_ref_to_sp_global = Qref * offset_ref_to_sp_local;

    std::array<Vec3, N> u;
    std::array<Vec3, N> theta;
    std::array<Mat3, N> R;
    std::array<Mat3, N> D;
    std::array<Vec3, N> x;
    std::array<std::array<Mat3, 3>, N> dR;

    for (Index node = 0; node < N; ++node) {
        u[node]     = q.template segment<3>(node * dofs_per_node);
        theta[node] = q.template segment<3>(node * dofs_per_node + 3);

        if (with_B) {
            math::so3::rotation_matrix_first_derivatives(
                theta[node], R[node], dR[node]);
        } else {
            R[node] = math::so3::rotation_matrix(theta[node]);
        }

        D[node] = R[node] * Qref;
        x[node] = node_position_reference(node)
                + u[node]
                + R[node] * offset_ref_to_sp_global;
    }

    const Vec3 centerline_derivative = (x[1] - x[0]) / L;

    // Geodesic Crisfield/Jelenic-style interpolation for the relative rotation.
    // D0^T D1 is invariant under a superposed spatial rigid-body rotation.
    const Mat3 relative_rotation = D[0].transpose() * D[1];
    const Vec3 relative_vector   = rotation_vector_from_matrix(relative_rotation);
    const Precision relative_angle = relative_vector.norm();

    logging::error(relative_angle < Precision(3.14159265358979323846) - Precision(1e-7),
        "B31: relative nodal rotation is too close to 180 degrees in element ",
        this->elem_id);

    const Vec3 curvature = relative_vector / L;
    const Mat3 Jrinv = with_B
        ? right_jacobian_inverse(relative_vector)
        : Mat3::Identity();

    Kinematics result;

    for (Index point = 0; point < N; ++point) {
        const Vec3 gamma = D[point].transpose() * centerline_derivative - e1_local;
        result.point[point].strain.template head<3>() = gamma;
        result.point[point].strain.template tail<3>() = curvature;

        if (!with_B) continue;

        for (Index column = 0; column < num_dofs; ++column) {
            const Index node = column / dofs_per_node;
            const Index dof  = column % dofs_per_node;
            const Precision sign = node == 0 ? Precision(-1) : Precision(1);

            Vec3 d_centerline = Vec3::Zero();
            Mat3 dD_point     = Mat3::Zero();
            Mat3 d_relative   = Mat3::Zero();

            if (dof < 3) {
                Vec3 direction = Vec3::Zero();
                direction(dof) = Precision(1);
                d_centerline = sign * direction / L;
            } else {
                const Index rotational_dof = dof - 3;
                const Mat3 dD = dR[node][rotational_dof] * Qref;
                const Vec3 dx = dR[node][rotational_dof] * offset_ref_to_sp_global;

                d_centerline = sign * dx / L;
                if (node == point) dD_point = dD;

                if (node == 0) {
                    d_relative = dD.transpose() * D[1];
                } else {
                    d_relative = D[0].transpose() * dD;
                }
            }

            const Vec3 d_gamma =
                dD_point.transpose() * centerline_derivative
                + D[point].transpose() * d_centerline;

            Vec3 d_curvature = Vec3::Zero();
            if (dof >= 3) {
                const Mat3 body_increment = relative_rotation.transpose() * d_relative;
                const Vec3 eta = vee_skew(body_increment);
                d_curvature = (Jrinv * eta) / L;
            }

            result.point[point].B.template block<3, 1>(0, column) = d_gamma;
            result.point[point].B.template block<3, 1>(3, column) = d_curvature;
        }
    }

    return result;
}

StaticMatrix<6, 6> B31::constitutive_matrix() {
    const auto* profile = get_profile();
    const auto* elastic = get_elasticity();

    const Precision EA   = elastic->youngs * profile->area_;
    const Precision GAy  = elastic->shear  * profile->shear_area_y_;
    const Precision GAz  = elastic->shear  * profile->shear_area_z_;
    const Precision GJ   = elastic->shear  * profile->torsion_inertia_;
    const Precision EIy  = elastic->youngs * profile->inertia_y_;
    const Precision EIz  = elastic->youngs * profile->inertia_z_;
    const Precision EIyz = elastic->youngs * profile->product_inertia_yz_;

    StaticMatrix<3, 3> A = StaticMatrix<3, 3>::Zero();
    A(0, 0) = EA;
    A(1, 1) = GAy;
    A(2, 2) = GAz;

    // The material/centroid point is SMP and the shear point is SP.
    // d = SMP - SP generates the Hodges extension/bending coupling.
    const Precision dy = -profile->offset_y_;
    const Precision dz = -profile->offset_z_;

    StaticMatrix<3, 3> S = StaticMatrix<3, 3>::Zero();
    S(0, 1) =  dz;
    S(0, 2) = -dy;

    StaticMatrix<3, 3> D = StaticMatrix<3, 3>::Zero();
    D(0, 0) = GJ;
    D(1, 1) = EIy;
    D(2, 2) = EIz;
    D(1, 2) = -EIyz;
    D(2, 1) = -EIyz;

    StaticMatrix<6, 6> C = StaticMatrix<6, 6>::Zero();
    C.template block<3, 3>(0, 0) = A;
    C.template block<3, 3>(0, 3) = A * S;
    C.template block<3, 3>(3, 0) = S.transpose() * A;
    C.template block<3, 3>(3, 3) = D + S.transpose() * A * S;
    return C;
}

Vec6 B31::thermal_strain(Index point, const Field* temperature) {
    Vec6 strain = Vec6::Zero();
    if (!temperature) return strain;

    auto material = get_material();
    if (!material->has_thermal_expansion()) return strain;

    const Precision value = (*temperature)(
        static_cast<Index>(node_ids[point]), 0);
    strain(0) =
        material->get_thermal_expansion()
        * (value - material->get_thermal_zero_temperature());
    return strain;
}

std::array<Vec6, B31::N> B31::section_resultants(
    const Kinematics& state,
    const Field* temperature
) {
    const StaticMatrix<6, 6> C = constitutive_matrix();
    std::array<Vec6, N> resultants;

    for (Index point = 0; point < N; ++point) {
        resultants[point] =
            C * (state.point[point].strain - thermal_strain(point, temperature));
    }
    return resultants;
}

B31::Vec12 B31::exact_internal_force(
    const Vec12& q,
    const Field* temperature
) {
    const Kinematics state = kinematics(q, true);
    const auto resultants = section_resultants(state, temperature);
    const Precision weight = reference_length() / Precision(2);

    Vec12 force = Vec12::Zero();
    for (Index point = 0; point < N; ++point) {
        force.noalias() +=
            weight * state.point[point].B.transpose() * resultants[point];
    }
    return force;
}

B31::Mat12 B31::material_tangent(const Kinematics& state) {
    const StaticMatrix<6, 6> C = constitutive_matrix();
    const Precision weight = reference_length() / Precision(2);

    Mat12 tangent = Mat12::Zero();
    for (Index point = 0; point < N; ++point) {
        tangent.noalias() +=
            weight * state.point[point].B.transpose() * C * state.point[point].B;
    }
    return tangent;
}

B31::Mat12 B31::geometric_tangent_from_resultants(
    const Vec12& q,
    const std::array<Vec6, N>& resultants
) {
    const Precision L     = reference_length();
    const Precision weight = L / Precision(2);
    const Mat3 Qref       = reference_frame();

    const Profile* profile = get_profile();
    const Vec3 offset_ref_to_sp_local(
        Precision(0),
        profile->offset_y_ - profile->reference_y_,
        profile->offset_z_ - profile->reference_z_);
    const Vec3 offset_ref_to_sp_global =
        Qref * offset_ref_to_sp_local;

    std::array<Vec3, N> theta;
    std::array<Mat3, N> R;
    std::array<Mat3, N> D;
    std::array<Vec3, N> x;
    std::array<std::array<Mat3, 3>, N> dR;
    std::array<std::array<std::array<Mat3, 3>, 3>, N> d2R;

    for (Index node = 0; node < N; ++node) {
        const Vec3 u = q.template segment<3>(node * dofs_per_node);
        theta[node] = q.template segment<3>(node * dofs_per_node + 3);

        math::so3::rotation_matrix_second_derivatives(
            theta[node],
            R[node],
            dR[node],
            d2R[node]
        );

        D[node] = R[node] * Qref;
        x[node] = node_position_reference(node)
                + u
                + R[node] * offset_ref_to_sp_global;
    }

    const Vec3 centerline_derivative = (x[1] - x[0]) / L;
    const Mat3 relative_rotation     = D[0].transpose() * D[1];
    const Vec3 relative_vector       =
        rotation_vector_from_matrix(relative_rotation);

    const Precision relative_angle = relative_vector.norm();
    logging::error(
        relative_angle < Precision(3.14159265358979323846) - Precision(1e-7),
        "B31: relative nodal rotation is too close to 180 degrees in element ",
        this->elem_id
    );

    const Mat3 Jrinv = right_jacobian_inverse(relative_vector);

    // First derivatives needed by the exact Hessian of the generalized strains.
    std::array<Vec3, num_dofs> d_centerline;
    std::array<Mat3, num_dofs> d_relative;
    std::array<Vec3, num_dofs> relative_increment;
    std::array<Vec3, num_dofs> d_relative_vector;

    for (Index column = 0; column < num_dofs; ++column) {
        d_centerline[column].setZero();
        d_relative[column].setZero();
        relative_increment[column].setZero();
        d_relative_vector[column].setZero();

        const Index node = column / dofs_per_node;
        const Index dof  = column % dofs_per_node;
        const Precision sign =
            node == 0 ? Precision(-1) : Precision(1);

        if (dof < 3) {
            Vec3 direction = Vec3::Zero();
            direction(dof) = Precision(1);
            d_centerline[column] = sign * direction / L;
            continue;
        }

        const Index a = dof - 3;
        const Mat3 dD = dR[node][a] * Qref;

        d_centerline[column] =
            sign * (dR[node][a] * offset_ref_to_sp_global) / L;

        if (node == 0) {
            d_relative[column] = dD.transpose() * D[1];
        } else {
            d_relative[column] = D[0].transpose() * dD;
        }

        relative_increment[column] =
            vee_skew(relative_rotation.transpose() * d_relative[column]);
        d_relative_vector[column] =
            Jrinv * relative_increment[column];
    }

    const auto first_director_derivative =
        [&](Index point, Index column) {
            Mat3 value = Mat3::Zero();

            const Index node = column / dofs_per_node;
            const Index dof  = column % dofs_per_node;

            if (dof >= 3 && node == point) {
                value = dR[node][dof - 3] * Qref;
            }
            return value;
        };

    const auto second_director_derivative =
        [&](Index point, Index column_a, Index column_b) {
            Mat3 value = Mat3::Zero();

            const Index node_a = column_a / dofs_per_node;
            const Index node_b = column_b / dofs_per_node;
            const Index dof_a  = column_a % dofs_per_node;
            const Index dof_b  = column_b % dofs_per_node;

            if (dof_a >= 3 && dof_b >= 3
                && node_a == point && node_b == point) {
                value = d2R[point][dof_a - 3][dof_b - 3] * Qref;
            }
            return value;
        };

    Mat12 geometric = Mat12::Zero();

    for (Index column_a = 0; column_a < num_dofs; ++column_a) {
        for (Index column_b = 0; column_b < num_dofs; ++column_b) {
            const Index node_a = column_a / dofs_per_node;
            const Index node_b = column_b / dofs_per_node;
            const Index dof_a  = column_a % dofs_per_node;
            const Index dof_b  = column_b % dofs_per_node;

            Vec3 d2_centerline = Vec3::Zero();
            Mat3 d2_relative   = Mat3::Zero();

            if (dof_a >= 3 && dof_b >= 3) {
                const Index a = dof_a - 3;
                const Index b = dof_b - 3;

                if (node_a == node_b) {
                    const Precision sign =
                        node_a == 0 ? Precision(-1) : Precision(1);

                    d2_centerline =
                        sign
                        * (d2R[node_a][a][b] * offset_ref_to_sp_global)
                        / L;

                    const Mat3 d2D = d2R[node_a][a][b] * Qref;
                    if (node_a == 0) {
                        d2_relative = d2D.transpose() * D[1];
                    } else {
                        d2_relative = D[0].transpose() * d2D;
                    }
                } else if (node_a == 0) {
                    d2_relative =
                        (dR[0][a] * Qref).transpose()
                        * (dR[1][b] * Qref);
                } else {
                    d2_relative =
                        (dR[0][b] * Qref).transpose()
                        * (dR[1][a] * Qref);
                }
            }

            // Exact second derivative of phi = Log(D0^T D1).
            const Vec3 d_increment =
                vee_skew(
                    d_relative[column_b].transpose() * d_relative[column_a]
                    + relative_rotation.transpose() * d2_relative
                );

            const Mat3 dJrinv =
                right_jacobian_inverse_directional(
                    relative_vector,
                    d_relative_vector[column_b]
                );

            const Vec3 d2_curvature =
                (dJrinv * relative_increment[column_a]
                 + Jrinv * d_increment) / L;

            for (Index point = 0; point < N; ++point) {
                const Mat3 dD_a =
                    first_director_derivative(point, column_a);
                const Mat3 dD_b =
                    first_director_derivative(point, column_b);
                const Mat3 d2D =
                    second_director_derivative(
                        point, column_a, column_b);

                const Vec3 d2_gamma =
                      d2D.transpose() * centerline_derivative
                    + dD_a.transpose() * d_centerline[column_b]
                    + dD_b.transpose() * d_centerline[column_a]
                    + D[point].transpose() * d2_centerline;

                Vec6 d2_strain = Vec6::Zero();
                d2_strain.template head<3>() = d2_gamma;
                d2_strain.template tail<3>() = d2_curvature;

                geometric(column_a, column_b) +=
                    weight * resultants[point].dot(d2_strain);
            }
        }
    }

    // The exact operator is a Hessian. Average only roundoff-level asymmetry
    // from the SO(3) logarithm/Jacobian evaluations.
    return Precision(0.5) * (geometric + geometric.transpose());
}

std::array<Vec6, B31::N> B31::continued_resultants(
    const Kinematics& base_state,
    const std::array<Vec6, N>& base_resultants,
    const std::array<Vec6, N>& target_temperature_resultants,
    const Vec12& delta
) {
    (void) base_resultants;

    const StaticMatrix<6, 6> C = constitutive_matrix();
    std::array<Vec6, N> resultants;

    for (Index point = 0; point < N; ++point) {
        resultants[point] =
            target_temperature_resultants[point]
            + C * (base_state.point[point].B * delta);
    }
    return resultants;
}

MapMatrix B31::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force_output,
    const Field* target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature,
    bool         update_state
) {
    (void) update_state;

    const bool with_tangent   = tangent_buffer != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force_output != nullptr;

    if (!with_tangent && !with_geometric && !with_force) {
        return MapMatrix(nullptr, 0, 0);
    }

    logging::error(!with_force || target_displacement != nullptr,
        "B31: internal force evaluation requires displacement");
    logging::error(!with_force || internal_force_output->components >= dofs_per_node,
        "B31: internal force requires six nodal components");

    const Vec12 q0 = element_displacement(base_displacement);
    const Vec12 q  = target_displacement
        ? element_displacement(target_displacement)
        : q0;
    const Vec12 delta = q - q0;

    const Kinematics base_state = kinematics(q0, true);
    const auto base_resultants =
        section_resultants(base_state, base_temperature);
    const auto target_temperature_resultants =
        section_resultants(base_state, target_temperature);

    const bool affine_force =
        with_force && target_displacement != base_displacement;
    const bool need_complete_tangent = with_tangent || affine_force;

    Mat12 complete = Mat12::Zero();
    if (need_complete_tangent) {
        complete = material_tangent(base_state)
                 + geometric_tangent_from_resultants(q0, base_resultants);
    }

    if (with_geometric) {
        const auto continued = continued_resultants(
            base_state,
            base_resultants,
            target_temperature_resultants,
            delta);

        std::array<Vec6, N> increment;
        for (Index point = 0; point < N; ++point) {
            increment[point] = continued[point] - base_resultants[point];
        }

        MapMatrix mapped(geometric_tangent_buffer, num_dofs, num_dofs);
        mapped = geometric_tangent_from_resultants(q0, increment);
    }

    if (with_force) {
        const Precision weight = reference_length() / Precision(2);
        Vec12 force = Vec12::Zero();

        for (Index point = 0; point < N; ++point) {
            force.noalias() +=
                weight
                * base_state.point[point].B.transpose()
                * target_temperature_resultants[point];
        }

        if (affine_force) {
            force.noalias() += complete * delta;
        }

        for (Index node = 0; node < N; ++node) {
            const Index node_id = static_cast<Index>(node_ids[node]);
            for (Index dof = 0; dof < dofs_per_node; ++dof) {
                (*internal_force_output)(node_id, dof) +=
                    force(node * dofs_per_node + dof);
            }
        }
    }

    if (with_tangent) {
        MapMatrix mapped(tangent_buffer, num_dofs, num_dofs);
        mapped = complete;
    }

    if (with_tangent) {
        return MapMatrix(tangent_buffer, num_dofs, num_dofs);
    }
    if (with_geometric) {
        return MapMatrix(geometric_tangent_buffer, num_dofs, num_dofs);
    }
    return MapMatrix(nullptr, 0, 0);
}

StaticMatrix<12, 12> B31::stiffness_impl() {
    const Vec12 q = Vec12::Zero();
    const Kinematics state = kinematics(q, true);
    return material_tangent(state);
}

StaticMatrix<12, 12> B31::stiffness_geom_impl(
    const Field& target_displacement,
    const Field* target_temperature,
    const Field* base_temperature
) {
    const Vec12 q0 = Vec12::Zero();
    const Vec12 q  = element_displacement(&target_displacement);
    const Kinematics base_state = kinematics(q0, true);

    const auto base_resultants =
        section_resultants(base_state, base_temperature);
    const auto target_temperature_resultants =
        section_resultants(base_state, target_temperature);
    const auto continued = continued_resultants(
        base_state,
        base_resultants,
        target_temperature_resultants,
        q);

    std::array<Vec6, N> increment;
    for (Index point = 0; point < N; ++point) {
        increment[point] = continued[point] - base_resultants[point];
    }

    return geometric_tangent_from_resultants(q0, increment);
}

StaticMatrix<12, 12> B31::mass_impl() {
    const auto material = get_material();
    if (!material->has_density()) {
        return Mat12::Zero();
    }

    const Precision rho = material->get_density();
    const Precision L   = reference_length();
    const Profile* profile = get_profile();
    const Mat3 Qref = reference_frame();

    StaticMatrix<3, 3> Jarea = StaticMatrix<3, 3>::Zero();
    Jarea(0, 0) = profile->inertia_y_ + profile->inertia_z_;
    Jarea(1, 1) = profile->inertia_y_;
    Jarea(2, 2) = profile->inertia_z_;
    Jarea(1, 2) = -profile->product_inertia_yz_;
    Jarea(2, 1) = -profile->product_inertia_yz_;

    StaticMatrix<6, 6> M_smp = StaticMatrix<6, 6>::Zero();
    M_smp.template block<3, 3>(0, 0) =
        rho * profile->area_ * Mat3::Identity();
    M_smp.template block<3, 3>(3, 3) = rho * Jarea;

    const Vec3 offset_ref_to_smp(
        Precision(0),
        -profile->reference_y_,
        -profile->reference_z_);

    StaticMatrix<6, 6> B = StaticMatrix<6, 6>::Identity();
    B.template block<3, 3>(0, 3) = -math::skew(offset_ref_to_smp);

    const StaticMatrix<6, 6> M_ref_local =
        B.transpose() * M_smp * B;

    StaticMatrix<6, 6> Tnode = StaticMatrix<6, 6>::Zero();
    Tnode.template block<3, 3>(0, 0) = Qref.transpose();
    Tnode.template block<3, 3>(3, 3) = Qref.transpose();

    const StaticMatrix<6, 6> M_ref_global =
        Tnode.transpose() * M_ref_local * Tnode;

    Mat12 mass = Mat12::Zero();
    const Precision nodal_weight = L / Precision(2);
    mass.template block<6, 6>(0, 0) = nodal_weight * M_ref_global;
    mass.template block<6, 6>(6, 6) = nodal_weight * M_ref_global;
    return mass;
}

RowMatrix B31::stress_strain_nodal_rst() {
    RowMatrix rst(N, 3);
    rst.setZero();
    rst(0, 0) = Precision(-1);
    rst(1, 0) = Precision( 1);
    return rst;
}

RowMatrix B31::stress_strain_ip_rst() {
    return stress_strain_nodal_rst();
}

Vec6 B31::resultant_about_reference(const Vec6& shear_point_resultant) const {
    const auto* profile = this->_section->template as<BeamSection>()->profile_.get();
    const Vec3 offset_ref_to_sp(
        Precision(0),
        profile->offset_y_ - profile->reference_y_,
        profile->offset_z_ - profile->reference_z_);

    Vec6 result = shear_point_resultant;
    result.template tail<3>() +=
        offset_ref_to_sp.cross(result.template head<3>());
    return result;
}

void B31::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     target_displacement,
    const Field*     target_temperature,
    const RowMatrix& rst,
    const Field*     base_displacement,
    const Field*     base_temperature
) {
    logging::error(strain != nullptr || stress != nullptr,
        "B31: compute_stress_strain requires at least one output field");
    logging::error((!strain || strain->domain == FieldDomain::ELEMENT_NODAL)
                && (!stress || stress->domain == FieldDomain::ELEMENT_NODAL),
        "B31: stress/strain recovery requires ELEMENT_NODAL output");

    const Vec12 q0 = element_displacement(base_displacement);
    const Vec12 q  = element_displacement(&target_displacement);
    const Vec12 delta = q - q0;

    const Kinematics base_state = kinematics(q0, true);
    const auto base_resultants =
        section_resultants(base_state, base_temperature);
    const auto target_temperature_resultants =
        section_resultants(base_state, target_temperature);
    const auto resultants = continued_resultants(
        base_state,
        base_resultants,
        target_temperature_resultants,
        delta);

    const Index offset = static_cast<Index>(this->elem_nodal_offset);

    for (Index row_local = 0; row_local < static_cast<Index>(rst.rows()); ++row_local) {
        const Index point = row_local == 0 ? 0 : 1;
        const Index row = offset + row_local;
        const Vec6 generalized_strain =
            base_state.point[point].strain
            + base_state.point[point].B * delta;
        const Vec6 generalized_resultant =
            resultant_about_reference(resultants[point]);

        if (strain) {
            for (Index component = 0; component < strain->components; ++component) {
                (*strain)(row, component) = Precision(0);
            }
            const Index count = std::min<Index>(strain->components, n_strain);
            for (Index component = 0; component < count; ++component) {
                (*strain)(row, component) = generalized_strain(component);
            }
        }

        if (stress) {
            for (Index component = 0; component < stress->components; ++component) {
                (*stress)(row, component) = Precision(0);
            }
            const Index count = std::min<Index>(stress->components, n_strain);
            for (Index component = 0; component < count; ++component) {
                (*stress)(row, component) = generalized_resultant(component);
            }
        }
    }
}

bool B31::compute_beam_section_forces(
    Field&       section_forces,
    const Field& target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    const Vec12 q0 = element_displacement(base_displacement);
    const Vec12 q  = element_displacement(&target_displacement);
    const Vec12 delta = q - q0;

    const Kinematics base_state = kinematics(q0, true);
    const auto base_resultants =
        section_resultants(base_state, base_temperature);
    const auto target_temperature_resultants =
        section_resultants(base_state, target_temperature);
    const auto resultants = continued_resultants(
        base_state,
        base_resultants,
        target_temperature_resultants,
        delta);

    const Index offset = static_cast<Index>(this->elem_nodal_offset);
    for (Index point = 0; point < N; ++point) {
        const Vec6 value = resultant_about_reference(resultants[point]);
        const Index row = offset + point;

        for (Index component = 0; component < section_forces.components; ++component) {
            section_forces(row, component) = Precision(0);
        }

        const Index count = std::min<Index>(section_forces.components, n_strain);
        for (Index component = 0; component < count; ++component) {
            section_forces(row, component) = value(component);
        }
    }

    return true;
}

} // namespace fem::model
