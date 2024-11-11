/******************************************************************************
 * @file element_solid.ipp
 * @brief Implementation of the SolidElement class template. This file contains
 * the definitions for methods declared in the SolidElement class, including
 * the computation of the strain-displacement matrix, Jacobian, stiffness and
 * mass matrices, and other element-level calculations.
 *
 * @date Created on 12.06.2023
 * @author Finn Eggers
 ******************************************************************************/

#pragma once

#include "element_solid.h"

namespace fem::model {

//-----------------------------------------------------------------------------
// compute_stress_strain_nodal
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_stress_strain_nodal(NodeData& displacement, NodeData& stress, NodeData& strain) {
    auto local_node_coords = this->node_coords_local();
    auto global_node_coords = this->node_coords_global();

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Dim n = 0; n < n_nodes(); n++) {
        Precision r = local_node_coords(n, 0);
        Precision s = local_node_coords(n, 1);
        Precision t = local_node_coords(n, 2);
        Precision det;
        auto node_id = node_ids[n];

        StaticMatrix<n_strain, D * N> B = this->strain_displacements(global_node_coords, r, s, t, det, false);

        // some elements may have inf or nan determinant. if so, extrapolate
        if (std::isnan(det) || std::isinf(det)) {
            det = 0;
        }
        logging::error(det < 1e10, "invalid determinant encountered in element ", elem_id,
                               "\ndet        : ", det,
                               "\nCoordinates: ", global_node_coords);

        if (det > 0) {
            StaticMatrix<n_strain, n_strain> E = material_matrix(r, s, t);

            auto strains = B * local_displacement;
            auto stresses = E * strains;

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
            NodeData ip_xyz{this->n_integration_points(), 3};
            NodeData ip_stress{this->n_integration_points(), 6};
            NodeData ip_strain{this->n_integration_points(), 6};
            ip_xyz.setZero();
            ip_stress.setZero();
            ip_strain.setZero();

            compute_stress_strain(displacement, ip_stress, ip_strain, ip_xyz);
            auto res1 =
                fem::math::interpolate::interpolate(ip_xyz,
                                                    ip_stress,
                                                    global_node_coords.row(n),
                                                    nullptr,
                                                    0.7,
                                                    fem::math::interpolate::InterpolationFunction::QUADRATIC);
            auto res2 =
                fem::math::interpolate::interpolate(ip_xyz,
                                                    ip_strain,
                                                    global_node_coords.row(n),
                                                    nullptr,
                                                    0.7,
                                                    fem::math::interpolate::InterpolationFunction::QUADRATIC);

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

//-----------------------------------------------------------------------------
// compute_stress_strain
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_stress_strain(NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) {
    auto global_node_coords = this->node_coords_global();

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    auto scheme = this->integration_scheme();

    for (Index n = 0; n < scheme.count(); n++) {
        Precision r = scheme.get_point(n).r;
        Precision s = scheme.get_point(n).s;
        Precision t = scheme.get_point(n).t;
        Precision det;

        StaticMatrix<N, 1> shape_func = this->shape_function(r, s, t);
        StaticMatrix<n_strain, D * N> B = this->strain_displacements(global_node_coords, r, s, t, det);
        StaticMatrix<n_strain, n_strain> E = material_matrix(r, s, t);

        auto strains = B * local_displacement;
        auto stresses = E * strains;

        Precision x = 0;
        Precision y = 0;
        Precision z = 0;

        for (Dim j = 0; j < n_strain; j++) {
            strain(n, j) = strains(j);
            stress(n, j) = stresses(j);
        }

        for (Index j = 0; j < N; j++) {
            x += shape_func(j) * global_node_coords(j, 0);
            y += shape_func(j) * global_node_coords(j, 1);
            z += shape_func(j) * global_node_coords(j, 2);
        }

        xyz(n, 0) = x;
        xyz(n, 1) = y;
        xyz(n, 2) = z;
    }
}

//-----------------------------------------------------------------------------
// compute_compliance
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_compliance(NodeData& displacement, ElementData& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(buffer);

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0) = strain_energy;
}

template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(NodeData& displacement, ElementData& result) {
    // the equation for the compliance is:
    // C = u^T * K * u
    // C = u^T * (B^T * E * B) * u
    // C = u^T * B^T * R^T * E * R * B * u
    // C = (B * u)^T * R^T * E * R * (B * u)
    // C = e^T * R^T * E * R * e        where e = B * u

    // Evaluating this must be done using integration over the element since K comes from an integration

    // compute the u-vector
    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    // check for angles
    if (!this->_model_data->elem_data.has(TOPO_ANGLES)) return;

    Vec3 angles = this->_model_data->elem_data.get(TOPO_ANGLES).row(this->elem_id);
    auto R = cos::RectangularSystem::euler(angles(0), angles(1), angles(2)).get_axes(Vec3(0,0,0));
    auto dR_d1 = cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2));
    auto dR_d2 = cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2));
    auto dR_d3 = cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2));

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

        // compute the E matrix
        StaticMatrix<n_strain, n_strain> E = material_matrix(r, s, t);

        // compute the strain vector
        StaticVector<n_strain> strain = B * local_displacement;

        // compute the stress vector
        StaticVector<n_strain> stress = E * strain;

        // derivative of the rotated material matrix (R^T E R) w.r.t each angle
        auto dCd1 = this->_material->elasticity()->get_transformed_derivative(R, dR_d1);
        auto dCd2 = this->_material->elasticity()->get_transformed_derivative(R, dR_d2);
        auto dCd3 = this->_material->elasticity()->get_transformed_derivative(R, dR_d3);

        derivative(0) += w * strain.dot(dCd1 * strain) * det;
        derivative(1) += w * strain.dot(dCd2 * strain) * det;
        derivative(2) += w * strain.dot(dCd3 * strain) * det;
    }
}

}  // namespace fem::model
