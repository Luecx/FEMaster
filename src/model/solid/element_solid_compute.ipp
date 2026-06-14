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

    auto node_coords = this->node_coords_global();
    auto local_displacement = this->nodal_data<3>(displacement);
    auto local_disp_mat = StaticMatrix<3, N>(local_displacement.transpose());
    auto local_displacement_vec = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Index n = 0; n < rst.rows(); ++n) {
        const Precision r = rst(n, 0);
        const Precision s = rst(n, 1);
        const Precision t = rst(n, 2);
        const Index row = static_cast<Index>(offset) + n;

        if (!use_green_lagrange_nl) {
            Precision det;
            StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(node_coords, r, s, t, det, false);
            if (det <= Precision(0) || std::isnan(det) || std::isinf(det)) {
                continue;
            }

            const Vec6 global_strain = B * local_displacement_vec;
            Vec6 local_strain;
            Vec6 local_stress;
            material_stress_strain(r, s, t, global_strain, local_stress, local_strain);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain) (*strain)(row, j) = local_strain(j);
                if (stress) (*stress)(row, j) = local_stress(j);
            }
            continue;
        }

        StaticMatrix<N, D> dN_local = this->shape_derivative(r, s, t);
        StaticMatrix<D, D> J0 = this->jacobian(node_coords, r, s, t);
        const Precision detJ0 = J0.determinant();
        logging::error(detJ0 > Precision(0),
                       "negative reference determinant encountered in nonlinear stress for element ", elem_id,
                       "\ndet        : ", detJ0,
                       "\nCoordinates: ", node_coords);

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

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain) (*strain)(row, j) = green_voigt(j);
            if (stress) (*stress)(row, j) = cauchy_voigt(j);
        }
    }
}

template<Index N>
void SolidElement<N>::compute_stress_state(Field& stress_state,
                                           const Field& displacement,
                                           int offset,
                                           bool use_green_lagrange_nl) {
    RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) {
        return;
    }
    compute_stress_strain(nullptr, &stress_state, displacement, rst, offset, use_green_lagrange_nl);
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

template<Index N>
bool SolidElement<N>::compute_shear_flow(Field& shear_flow,
                                         const Field& displacement,
                                         int offset) {
    (void) shear_flow;
    (void) displacement;
    (void) offset;
    return false;
}

template<Index N>
bool SolidElement<N>::compute_beam_section_forces(Field& section_forces,
                                                  const Field& displacement,
                                                  int offset) {
    (void) section_forces;
    (void) displacement;
    (void) offset;
    return false;
}

template<Index N>
bool SolidElement<N>::compute_shell_section_forces(Field& section_forces,
                                                   Field& contribution_count,
                                                   const Field& displacement) {
    (void) section_forces;
    (void) contribution_count;
    (void) displacement;
    return false;
}

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
