/**
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 */

#pragma once

#include "../../cos/rectangular_system.h"

namespace fem::model {

namespace detail_solid_nonlinear {

inline Mat3 voigt_to_tensor(const Vec6& v) {
    Mat3 tensor;
    tensor << v(0), v(5), v(4),
              v(5), v(1), v(3),
              v(4), v(3), v(2);
    return tensor;
}

inline Vec6 tensor_to_voigt(const Mat3& tensor) {
    Vec6 v;
    v << tensor(0, 0),
         tensor(1, 1),
         tensor(2, 2),
         tensor(1, 2),
         tensor(2, 0),
         tensor(0, 1);
    return v;
}

inline Vec6 green_lagrange_to_voigt(const Mat3& strain) {
    Vec6 v;
    v << strain(0, 0),
         strain(1, 1),
         strain(2, 2),
         Precision(2) * strain(1, 2),
         Precision(2) * strain(2, 0),
         Precision(2) * strain(0, 1);
    return v;
}

} // namespace detail_solid_nonlinear


template<Index N>
Strains SolidElement<N>::strain(Field& displacement, std::vector<Vec3>& rst) {
    auto global_node_coords = this->node_coords_global();
    auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    return {};
}

template<Index N>
Stresses SolidElement<N>::stress(Field& displacement, std::vector<Vec3>& rst) {
    return {};
}


//-----------------------------------------------------------------------------
// compute_stress_strain_nodal
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_stress_strain_nodal(Field& displacement, Field& stress, Field& strain) {
    auto local_node_coords = this->node_coords_local();
    auto global_node_coords = this->node_coords_global();

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Dim n = 0; n < n_nodes(); n++) {
        Precision r = local_node_coords(n, 0);
        Precision s = local_node_coords(n, 1);
        Precision t = local_node_coords(n, 2);
        Precision det;
        const Index node_id = static_cast<Index>(node_ids[n]);

        StaticMatrix<n_strain, D * N> B = this->strain_displacements(global_node_coords, r, s, t, det, false);

        // some elements may have inf or nan determinant. if so, extrapolate
        if (std::isnan(det) || std::isinf(det)) {
            det = 0;
        }
        logging::error(det < 1e10, "invalid determinant encountered in element ", elem_id,
                               "\ndet        : ", det,
                               "\nCoordinates: ", global_node_coords);

        if (det > 0) {
            auto global_strains = B * local_displacement;
            StaticVector<n_strain> strains;
            StaticVector<n_strain> stresses;
            material_stress_strain(r, s, t, global_strains, stresses, strains);

            // if any nan
            if (stresses.hasNaN()) {
                logging::error(false, "invalid stress encountered in element ", elem_id,
                               "\nStrains: ", strains,
                               "\nStresses: ", stresses);
            }

            for (int j = 0; j < n_strain; j++) {
                // check for any inf
                if (std::isinf(strains(j)) || std::isinf(stresses(j))) {
                        logging::error(false, "invalid stress encountered in element ", elem_id,
                   "\nStrains: ", strains,
                   "\nStresses: ", stresses);
                }

                strain(node_id, j) += strains(j);
                stress(node_id, j) += stresses(j);
            }
        } else {
            const Index n_ip = static_cast<Index>(this->n_integration_points());
            RowMatrix ip_xyz(n_ip, 3);
            ip_xyz.setZero();
            Field ip_stress{"IP_STRESS_LOCAL", FieldDomain::ELEMENT_IP, n_ip, 6};
            Field ip_strain{"IP_STRAIN_LOCAL", FieldDomain::ELEMENT_IP, n_ip, 6};
            ip_stress.set_zero();
            ip_strain.set_zero();

            // fill ip_xyz
            auto scheme = this->integration_scheme();
            for (int i = 0; i < scheme.count(); i++) {
                Precision r_ip = scheme.get_point(i).r;
                Precision s_ip = scheme.get_point(i).s;
                Precision t_ip = scheme.get_point(i).t;
                Precision w_ip = scheme.get_point(i).w;

                Precision x = 0;
                Precision y = 0;
                Precision z = 0;

                auto shape_func = this->shape_function(r_ip, s_ip, t_ip);
                for (Index j = 0; j < N; j++) {
                    x += shape_func(j) * global_node_coords(j, 0);
                    y += shape_func(j) * global_node_coords(j, 1);
                    z += shape_func(j) * global_node_coords(j, 2);
                }
                ip_xyz(i, 0) = x;
                ip_xyz(i, 1) = y;
                ip_xyz(i, 2) = z;
            }

            compute_stress_strain(ip_stress, ip_strain, displacement, 0);

            RowMatrix ip_stress_mat =
                Eigen::Map<const RowMatrix>(ip_stress.data(), ip_stress.rows, ip_stress.components);
            RowMatrix ip_strain_mat =
                Eigen::Map<const RowMatrix>(ip_strain.data(), ip_strain.rows, ip_strain.components);
            auto res1 =
                fem::math::interpolate::interpolate(ip_xyz,
                                                    ip_stress_mat,
                                                    global_node_coords.row(n),
                                                    nullptr,
                                                    0.7,
                                                    fem::math::interpolate::InterpolationFunction::CONSTANT);
            auto res2 =
                fem::math::interpolate::interpolate(ip_xyz,
                                                    ip_strain_mat,
                                                    global_node_coords.row(n),
                                                    nullptr,
                                                    0.7,
                                                    fem::math::interpolate::InterpolationFunction::CONSTANT);

            for (int j = 0; j < n_strain; j++) {

                // if any nan or inf values produced, error
                if (std::isnan(res1(j)) || std::isinf(res1(j))) {
                    logging::error(false, "invalid stress encountered in element ", elem_id,
                                   "\nStrains: ", res2,
                                   "\nStresses: ", res1);
                }

                stress(node_id, j) += res1(j);
                strain(node_id, j) += res2(j);
            }
        }
    }
}

template<Index N>
void SolidElement<N>::compute_stress_strain(Field& ip_stress, Field& ip_strain, Field& displacement, int ip_offset) {
    auto global_node_coords = this->node_coords_global();
    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
    auto scheme = this->integration_scheme();
    for (Index n = 0; n < scheme.count(); n++) {
        Precision r = scheme.get_point(n).r;
        Precision s = scheme.get_point(n).s;
        Precision t = scheme.get_point(n).t;
        Precision det;
        StaticMatrix<N, 1> shape_func       = this->shape_function(r, s, t);
        StaticMatrix<n_strain, D * N> B     = this->strain_displacements(global_node_coords, r, s, t, det);
        auto global_strains = B * local_displacement;
        StaticVector<n_strain> strains;
        StaticVector<n_strain> stresses;
        material_stress_strain(r, s, t, global_strains, stresses, strains);
        for (Dim j = 0; j < n_strain; j++) {
            const Index row = static_cast<Index>(ip_offset) + n;
            ip_stress(row, j) = stresses(j);
            ip_strain(row, j) = strains(j);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_ip_stress_nonlinear(Field& ip_stress, Field& displacement, int ip_offset) {
    auto reference_coords = this->node_coords_global();
    auto local_displacement = this->nodal_data<3>(displacement);
    auto scheme = this->integration_scheme();

    for (Index n = 0; n < scheme.count(); n++) {
        const Precision r = scheme.get_point(n).r;
        const Precision s = scheme.get_point(n).s;
        const Precision t = scheme.get_point(n).t;
        const Index row = static_cast<Index>(ip_offset) + n;

        StaticMatrix<N, D> dN_local = this->shape_derivative(r, s, t);
        StaticMatrix<D, D> J0 = this->jacobian(reference_coords, r, s, t);
        const Precision detJ0 = J0.determinant();
        logging::error(detJ0 > Precision(0),
                       "negative reference determinant encountered in nonlinear stress for element ", elem_id,
                       "\ndet        : ", detJ0,
                       "\nCoordinates: ", reference_coords);

        StaticMatrix<N, D> dN_dX = (J0.inverse() * dN_local.transpose()).transpose();

        Mat3 grad_u = Mat3::Zero();
        for (Index a = 0; a < N; ++a) {
            for (Index i = 0; i < D; ++i) {
                for (Index j = 0; j < D; ++j) {
                    grad_u(i, j) += local_displacement(a, i) * dN_dX(a, j);
                }
            }
        }

        const Mat3 F = Mat3::Identity() + grad_u;
        const Precision J = F.determinant();
        logging::error(J > Precision(0),
                       "non-positive deformation gradient determinant in element ", elem_id,
                       "\nJ: ", J,
                       "\nF: ", F);

        const Mat3 green_lagrange =
            Precision(0.5) * (F.transpose() * F - Mat3::Identity());
        const Vec6 green_voigt =
            detail_solid_nonlinear::green_lagrange_to_voigt(green_lagrange);

        const Vec6 second_pk_voigt = material_matrix(r, s, t) * green_voigt;

        const Mat3 second_pk = detail_solid_nonlinear::voigt_to_tensor(second_pk_voigt);
        const Mat3 cauchy = (F * second_pk * F.transpose()) / J;
        const Vec6 cauchy_voigt = detail_solid_nonlinear::tensor_to_voigt(cauchy);

        for (Dim j = 0; j < n_strain; j++) {
            ip_stress(row, j) = cauchy_voigt(j);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_internal_force_nonlinear(Field& node_forces,
                                                       const Field& ip_stress,
                                                       int ip_offset) {
    auto current_coords = this->node_coords_global();
    auto scheme = this->integration_scheme();

    for (Index n = 0; n < scheme.count(); n++) {
        const Precision r = scheme.get_point(n).r;
        const Precision s = scheme.get_point(n).s;
        const Precision t = scheme.get_point(n).t;
        const Precision w = scheme.get_point(n).w;
        const Index row = static_cast<Index>(ip_offset) + n;

        Precision det;
        StaticMatrix<n_strain, D * N> B =
            this->strain_displacements(current_coords, r, s, t, det);
        const Vec6 sigma = ip_stress.row_vec6(row);
        const StaticVector<D * N> local_force = B.transpose() * sigma * (det * w);

        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Index d = 0; d < D; ++d) {
                node_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }
}


// //-----------------------------------------------------------------------------
// // compute_stress_strain
// //-----------------------------------------------------------------------------
// template<Index N>
// void
// SolidElement<N>::compute_stress_strain(Field& displacement, Field& stress, Field& strain, Field& xyz) {
//     auto global_node_coords = this->node_coords_global();
//
//     auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
//     auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
//
//     auto scheme = this->integration_scheme();
//
//     for (Index n = 0; n < scheme.count(); n++) {
//         Precision r = scheme.get_point(n).r;
//         Precision s = scheme.get_point(n).s;
//         Precision t = scheme.get_point(n).t;
//         Precision det;
//
//         StaticMatrix<N, 1> shape_func = this->shape_function(r, s, t);
//         StaticMatrix<n_strain, D * N> B = this->strain_displacements(global_node_coords, r, s, t, det);
//         StaticMatrix<n_strain, n_strain> E = material_matrix(r, s, t);
//
//         auto strains = B * local_displacement;
//         auto stresses = E * strains;
//
//         Precision x = 0;
//         Precision y = 0;
//         Precision z = 0;
//
//         for (Dim j = 0; j < n_strain; j++) {
//             strain(n, j) = strains(j);
//             stress(n, j) = stresses(j);
//         }
//
//         for (Index j = 0; j < N; j++) {
//             x += shape_func(j) * global_node_coords(j, 0);
//             y += shape_func(j) * global_node_coords(j, 1);
//             z += shape_func(j) * global_node_coords(j, 2);
//         }
//
//         xyz(n, 0) = x;
//         xyz(n, 1) = y;
//         xyz(n, 2) = z;
//     }
// }

//-----------------------------------------------------------------------------
// compute_compliance
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(buffer);

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0) = strain_energy;
}

template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    // the equation for the compliance is:
    // C = u^T * K * u
    // C = u^T * (B^T * E * B) * u
    // C = u^T * B^T * R^T * E * R * B * u
    // C = (B * u)^T * R^T * E * R * (B * u)
    // C = e^T * R^T * E * R * e        where e = B * u
    // C = e^T (D) e
    // where D = R^T * E * R
    // it follows
    // dC/dtheta = de^T/dtheta * D * e + e^T * D * de/dtheta

    // Evaluating this must be done using integration over the element since K comes from an integration

    // compute the u-vector
    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    // check for angles
    if (!this->_model_data || !this->_model_data->material_orientation) return;

    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
                   "Field '", angles_field->name, "': material orientation requires 3 components");
    const Index row = static_cast<Index>(this->elem_id);
    Vec3 angles = angles_field->row_vec3(row);
    auto R_topo = cos::RectangularSystem::euler(angles(0), angles(1), angles(2)).get_axes(Vec3(0,0,0));
    auto dR_topo_d1 = cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2));
    auto dR_topo_d2 = cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2));
    auto dR_topo_d3 = cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2));

    Vec3 derivative = Vec3::Zero();

    // go through each integration point
    for (Index n = 0; n < this->integration_scheme().count(); n++) {
        Precision r = this->integration_scheme().get_point(n).r;
        Precision s = this->integration_scheme().get_point(n).s;
        Precision t = this->integration_scheme().get_point(n).t;
        Precision w = this->integration_scheme().get_point(n).w;
        Precision det;

        // compute the B matrix
        StaticMatrix<n_strain, D * N> B = this->strain_displacements(this->node_coords_global(), r, s, t, det);

        // compute the strain vector
        StaticVector<n_strain> strain = B * local_displacement;

        const Mat3 R_section = section_orientation_basis(r, s, t);
        const Mat3 R = R_section * R_topo;
        const Mat3 dR_d1 = R_section * dR_topo_d1;
        const Mat3 dR_d2 = R_section * dR_topo_d2;
        const Mat3 dR_d3 = R_section * dR_topo_d3;

        // derivative of the rotated material matrix with section orientation fixed
        auto dCd1 = this->material()->elasticity()->template get_transformed_derivative<3>(R, dR_d1);
        auto dCd2 = this->material()->elasticity()->template get_transformed_derivative<3>(R, dR_d2);
        auto dCd3 = this->material()->elasticity()->template get_transformed_derivative<3>(R, dR_d3);

        derivative(0) += w * strain.dot(dCd1 * strain) * det;
        derivative(1) += w * strain.dot(dCd2 * strain) * det;
        derivative(2) += w * strain.dot(dCd3 * strain) * det;
    }
    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

}  // namespace fem::model
