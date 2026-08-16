/**
 * @file element_solid_compute.ipp
 * @brief Implements solid recovery and nonlinear constitutive assembly.
 *
 * @author Finn Eggers
 * @date 16.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"

#include <cmath>

namespace fem::model {

template<Index N>
void SolidElement<N>::compute_stress_strain(Field* strain,
                                            Field* stress,
                                            const Field& displacement,
                                            const RowMatrix& rst,
                                            int offset,
                                            bool use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
        "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "SolidElement: stress/strain evaluation coordinates require at least 3 columns");

    const auto reference_coords = this->node_coords_reference();
    const auto current_coords = this->node_coords_current();
    const auto local_displacement = this->nodal_data<3>(displacement);
    const auto local_disp_mat = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_displacement_vec =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
    const auto& state_scheme = this->integration_scheme_stiffness();

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r = rst(n, 0);
        const Precision s = rst(n, 1);
        const Precision t = rst(n, 2);
        const Index row = static_cast<Index>(offset + n);

        Index state_ip = 0;
        auto point0 = state_scheme.get_point(0);
        Precision best = (r - point0.r) * (r - point0.r)
                       + (s - point0.s) * (s - point0.s)
                       + (t - point0.t) * (t - point0.t);
        for (Index ip = 1; ip < state_scheme.count(); ++ip) {
            const auto point = state_scheme.get_point(ip);
            const Precision distance = (r - point.r) * (r - point.r)
                                     + (s - point.s) * (s - point.s)
                                     + (t - point.t) * (t - point.t);
            if (distance < best) {
                best = distance;
                state_ip = ip;
            }
        }

        const Precision* committed_state = material_state_old(state_ip);

        if (!use_green_lagrange_nl) {
            Precision det0;
            const auto B = this->strain_displacements(
                reference_coords, r, s, t, det0, false);
            if (det0 <= Precision(0) || !std::isfinite(det0)) continue;

            const VolumeStrainLinearized global_strain(
                Vec6(B * local_displacement_vec));
            VolumeStressCauchy global_stress;
            recover_material(r, s, t, global_strain, committed_state, global_stress);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain) (*strain)(row, j) = global_strain.voigt()(j);
                if (stress) (*stress)(row, j) = global_stress.voigt()(j);
            }
            continue;
        }

        const Mat3 F = this->deformation_gradient(
            reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange global_strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        recover_material(r, s, t, global_strain, committed_state, second_pk);
        const VolumeStressCauchy cauchy = second_pk.to_cauchy(F);

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain) (*strain)(row, j) = global_strain.voigt()(j);
            if (stress) (*stress)(row, j) = cauchy.voigt()(j);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_stress_state(Field& stress_state,
                                           const Field& displacement,
                                           int offset,
                                           bool use_green_lagrange_nl) {
    const RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) return;

    const auto reference_coords = this->node_coords_reference();
    const auto current_coords = this->node_coords_current();
    const auto local_displacement = this->nodal_data<3>(displacement);
    const auto local_disp_mat = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_displacement_vec =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r = rst(n, 0);
        const Precision s = rst(n, 1);
        const Precision t = rst(n, 2);
        const Index row = static_cast<Index>(offset + n);
        const Index ip = static_cast<Index>(n);
        const Precision* state_old = material_state_old(ip);
        Precision* state_new = material_state_new(ip);

        if (!use_green_lagrange_nl) {
            Precision det0;
            const auto B = this->strain_displacements(
                reference_coords, r, s, t, det0, false);
            if (det0 <= Precision(0) || !std::isfinite(det0)) continue;

            const VolumeStrainLinearized eps(Vec6(B * local_displacement_vec));
            VolumeStressCauchy sigma;
            Mat6 tangent;
            update_material(r, s, t, eps, state_old, state_new, sigma, tangent);
            for (Dim component = 0; component < n_strain; ++component) {
                stress_state(row, component) = sigma.voigt()(component);
            }
            continue;
        }

        const Mat3 F = this->deformation_gradient(
            reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange eps =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        Mat6 tangent;
        update_material(r, s, t, eps, state_old, state_new, second_pk, tangent);
        for (Dim component = 0; component < n_strain; ++component) {
            stress_state(row, component) = second_pk.voigt()(component);
        }
    }
}

template<Index N>
MapMatrix SolidElement<N>::stiffness_tangent(Precision* buffer,
                                             Field& ip_stress_state,
                                             NodeData& nodal_forces,
                                             const Field& displacement) {
    logging::error(ip_stress_state.components >= n_strain,
        "SolidElement: nonlinear stress state requires at least six components");
    logging::error(nodal_forces.components >= D,
        "SolidElement: nonlinear internal force requires at least three nodal components");

    const auto reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();
    const auto local_displacement = this->nodal_data<3>(displacement);
    const auto local_disp_mat = StaticMatrix<3, N>(local_displacement.transpose());
    const auto u = Eigen::Map<const StaticVector<3 * N>>(
        local_disp_mat.data(), 3 * N);

    StaticMatrix<D * N, D * N> tangent =
        StaticMatrix<D * N, D * N>::Zero();

    if (!this->_model_data->geometric_nonlinearity) {
        // Material-nonlinear, geometrically linear assembly. J2 is integrated
        // here with infinitesimal strain, Cauchy stress and no geometric tangent.
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            Precision det0;
            const auto dN_dX = this->shape_derivatives_reference(
                reference_coords, point.r, point.s, point.t, det0);
            const auto B = this->strain_displacement(dN_dX);
            const VolumeStrainLinearized eps(Vec6(B * u));

            const Precision* state_old = material_state_old(ip);
            Precision* state_new = material_state_new(ip);
            VolumeStressCauchy sigma;
            Mat6 C_alg;
            update_material(
                point.r, point.s, point.t,
                eps, state_old, state_new, sigma, C_alg);

            const Index stress_row = this->ip_index(ip);
            for (Index component = 0; component < n_strain; ++component) {
                ip_stress_state(stress_row, component) = sigma.voigt()(component);
            }

            const Precision measure = point.w * det0;
            tangent.noalias() += measure * B.transpose() * C_alg * B;
            const StaticVector<D * N> local_force =
                B.transpose() * sigma.voigt() * measure;
            for (Index a = 0; a < N; ++a) {
                const Index node_id = static_cast<Index>(node_ids[a]);
                for (Dim d = 0; d < D; ++d) {
                    nodal_forces(node_id, d) += local_force(D * a + d);
                }
            }
        }

        tangent = Precision(0.5) * (tangent + tangent.transpose());
        MapMatrix mapped{buffer, D * N, D * N};
        mapped = tangent;
        return mapped;
    }

    // Existing Total-Lagrangian path for purely elastic finite-strain laws.
    const auto current_coords = this->node_coords_current();
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        Precision det0;
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Mat3 F = this->deformation_gradient(
            reference_coords, current_coords, point.r, point.s, point.t);
        const auto B = this->green_lagrange_strain_displacement(dN_dX, F);
        const VolumeStrainGreenLagrange eps =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        const Precision* state_old = material_state_old(ip);
        Precision* state_new = material_state_new(ip);
        VolumeStressPK2 stress;
        Mat6 C;
        update_material(
            point.r, point.s, point.t,
            eps, state_old, state_new, stress, C);

        const Index stress_row = this->ip_index(ip);
        for (Index component = 0; component < n_strain; ++component) {
            ip_stress_state(stress_row, component) = stress.voigt()(component);
        }

        const Precision measure = point.w * det0;
        tangent.noalias() += measure * B.transpose() * C * B;

        const Mat3 S = stress.tensor();
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();
            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision value = dNa.dot(S * dNb) * measure;
                for (Dim d = 0; d < D; ++d) {
                    tangent(D * a + d, D * b + d) += value;
                }
            }
        }

        const StaticVector<D * N> local_force =
            B.transpose() * stress.voigt() * measure;
        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Dim d = 0; d < D; ++d) {
                nodal_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }

    tangent = Precision(0.5) * (tangent + tangent.transpose());
    MapMatrix mapped{buffer, D * N, D * N};
    mapped = tangent;
    return mapped;
}

template<Index N>
void SolidElement<N>::compute_internal_force_nonlinear(Field& node_forces,
                                                       const Field& ip_stress) {
    const auto reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();

    if (!this->_model_data->geometric_nonlinearity) {
        for (Index ip = 0; ip < scheme.count(); ++ip) {
            const auto point = scheme.get_point(ip);
            Precision det0;
            const auto dN_dX = this->shape_derivatives_reference(
                reference_coords, point.r, point.s, point.t, det0);
            const auto B = this->strain_displacement(dN_dX);
            const auto sigma = ip_stress.row_vec6(this->ip_index(ip));
            const auto local_force =
                B.transpose() * sigma * (point.w * det0);
            for (Index a = 0; a < N; ++a) {
                const Index node_id = static_cast<Index>(node_ids[a]);
                for (Dim d = 0; d < D; ++d) {
                    node_forces(node_id, d) += local_force(D * a + d);
                }
            }
        }
        return;
    }

    const auto current_coords = this->node_coords_current();
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        Precision det0;
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto F = this->deformation_gradient(
            reference_coords, current_coords, point.r, point.s, point.t);
        const auto B = this->green_lagrange_strain_displacement(dN_dX, F);
        const auto S = ip_stress.row_vec6(this->ip_index(ip));
        const auto local_force = B.transpose() * S * (point.w * det0);
        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Dim d = 0; d < D; ++d) {
                node_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }
}

template<Index N>
void SolidElement<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(buffer);
    auto local_disp_mat = StaticMatrix<3, N>(
        this->nodal_data<3>(displacement).transpose());
    auto u = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
    result(elem_id, 0) = u.dot(K * u);
}

template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement,
                                                          Field& result) {
    if (!this->_model_data || !this->_model_data->material_orientation) return;
    logging::error(!this->_section->material_->has_plasticity(),
        "Compliance angle derivative is not defined for plastic material response");

    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name,
        "': material orientation requires 3 components");

    const Vec3 angles = angles_field->row_vec3(static_cast<Index>(this->elem_id));
    const Mat3 rotation = cos::RectangularSystem::euler(
        angles(0), angles(1), angles(2)).get_axes(Vec3::Zero());
    const std::array<Mat3, 3> rotation_derivatives {
        cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2))
    };

    auto local_disp_mat = StaticMatrix<3, N>(
        this->nodal_data<3>(displacement).transpose());
    const auto u = Eigen::Map<const StaticVector<3 * N>>(
        local_disp_mat.data(), 3 * N);
    const auto reference_coords = this->node_coords_reference();
    const auto& scheme = this->integration_scheme_stiffness();
    const Precision scaling = element_stiffness_scale();
    Vec3 derivative = Vec3::Zero();

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        Precision det0;
        const auto dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const auto B = this->strain_displacement(dN_dX);
        const StaticVector<n_strain> eps = B * u;
        const Vec3 position_reference = this->interpolate<D>(
            reference_coords, point.r, point.s, point.t);
        const auto derivatives = get_section()->tangent_rotation_derivatives(
            position_reference, rotation, rotation_derivatives);
        for (Index i = 0; i < 3; ++i) {
            derivative(i) += scaling * point.w * det0
                           * eps.dot(derivatives[i] * eps);
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

} // namespace fem::model
