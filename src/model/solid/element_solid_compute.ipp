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


}  // namespace fem::model
