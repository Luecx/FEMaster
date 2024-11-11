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

#include "../../cos/rectangular_system.h"

namespace fem::model {

//-----------------------------------------------------------------------------
// Interpolation
//-----------------------------------------------------------------------------
template<Index N>
template<Dim K>
StaticVector<K> SolidElement<N>::interpolate(StaticMatrix<N, K> values, Precision r, Precision s, Precision t) {
    StaticVector<N> wgt = shape_function(r,s,t);
    StaticVector<K> res {};
    res.setZero();
    for (Index i = 0; i < K; i++) {
        res += wgt(i) * values.col(i);
    }
    return res;
}

//-----------------------------------------------------------------------------
// strain_displacement
//-----------------------------------------------------------------------------
template<Index N>
StaticMatrix<SolidElement<N>::n_strain, SolidElement<N>::D * N>
SolidElement<N>::strain_displacement(const StaticMatrix<N, D>& shape_der_global) {
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();
    
    for (Index j = 0; j < N; j++) {
        Dim r1 = j * 3;
        Dim r2 = r1 + 1;
        Dim r3 = r1 + 2;
        
        B(0, r1) = shape_der_global(j, 0);  // dx/dr
        B(1, r2) = shape_der_global(j, 1);  // dy/ds
        B(2, r3) = shape_der_global(j, 2);  // dz/dt

        B(3, r2) = shape_der_global(j, 2);  // dy/dr
        B(3, r3) = shape_der_global(j, 1);  // dz/ds

        B(4, r1) = shape_der_global(j, 2);  // dz/dr
        B(4, r3) = shape_der_global(j, 0);  // dx/dt

        B(5, r1) = shape_der_global(j, 1);  // dx/ds
        B(5, r2) = shape_der_global(j, 0);  // dy/dt
    }
    
    return B;
}

//-----------------------------------------------------------------------------
// strain_displacements
//-----------------------------------------------------------------------------
template<Index N>
StaticMatrix<SolidElement<N>::n_strain, SolidElement<N>::D * N>
SolidElement<N>::strain_displacements(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t, Precision& det, bool check_det) {
    StaticMatrix<N, D> local_shape_der = shape_derivative(r, s, t);
    StaticMatrix<D, D> jac = jacobian(node_coords, r, s, t);
    
    det = jac.determinant();
    StaticMatrix<D, D> inv = jac.inverse();
    
    if (check_det) {
        logging::error(det > 0, "negative determinant encountered in element ", elem_id,
                       "\ndet        : ", det,
                       "\nCoordinates: ", node_coords,
                       "\nJacobi     : ", jac);
    }
    
    StaticMatrix<N, D> global_shape_der = (inv * local_shape_der.transpose()).transpose();
    return strain_displacement(global_shape_der);
}

//-----------------------------------------------------------------------------
// jacobian
//-----------------------------------------------------------------------------
template<Index N>
StaticMatrix<SolidElement<N>::D, SolidElement<N>::D>
SolidElement<N>::jacobian(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t) {
    StaticMatrix<N, D> local_shape_derivative = shape_derivative(r, s, t);
    StaticMatrix<D, D> jacobian {};
    
    for (Dim m = 0; m < D; m++) {
        for (Dim n = 0; n < D; n++) {
            Precision dxn_drm = 0;
            for (Dim k = 0; k < N; k++) {
                dxn_drm += node_coords(k, n) * local_shape_derivative(k, m);
            }
            jacobian(m, n) = dxn_drm;
        }
    }
    
    return jacobian;
}

template<Index N>
StaticMatrix<SolidElement<N>::n_strain, SolidElement<N>::n_strain>
SolidElement<N>::mat_matrix(const StaticMatrix<N, D>& node_coords, Precision r, Precision s, Precision t) {
    logging::error(_elem_data_dict != nullptr, "no _elem_data_dict assigned to element ", elem_id);
    logging::error(_node_data_dict != nullptr, "no _node_data_dict assigned to element ", elem_id);
    logging::error(_material != nullptr, "no _material assigned to element ", elem_id);
    logging::error(_material->has_elasticity(), "_material has no elasticity components assigned at element ", elem_id);

    // get the orientation
    if (_elem_data_dict->has(MAT_ANGLES)) {
        Vec3 orient_angles = _elem_data_dict->get(MAT_ANGLES).row(this->elem_id);
        Vec3 location      = this->interpolate(node_coords);
        cos::RectangularSystem local_system(orient_angles(0), orient_angles(1), orient_angles(2));
        return _material->elasticity()->template get_transformed<D>(local_system.get_axes(Vec3(0,0,0)));
    }
    return _material->elasticity()->template get<D>();
}


//-----------------------------------------------------------------------------
// nodal_data
//-----------------------------------------------------------------------------
template<Index N>
template<Dim K>
StaticMatrix<N, K>
SolidElement<N>::nodal_data(const NodeData& full_data, Index offset, Index stride) {
    StaticMatrix<N, K> res {};
    runtime_assert(full_data.cols() >= offset + stride * D, "cannot extract this many elements from the data");
    
    for (Dim m = 0; m < N; m++) {
        for (Dim j = 0; j < K; j++) {
            Index n = j * stride + offset;
            res(m, j) = full_data(node_ids[m], n);
        }
    }
    
    return res;
}

//-----------------------------------------------------------------------------
// dofs
//-----------------------------------------------------------------------------
template<Index N>
ElDofs
SolidElement<N>::dofs() {
    return ElDofs {true, true, true, false, false, false};
}

//-----------------------------------------------------------------------------
// stiffness
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::stiffness(NodeData& position, Precision* buffer) {
    logging::error(_material != nullptr, "no _material assigned to element ", elem_id);
    logging::error(_material->has_elasticity(), "_material has no elasticity components assigned at element ", elem_id);

    StaticMatrix<N, D> node_coords = this->node_coords_global(position);

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &node_coords](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det;
            StaticMatrix<n_strain, n_strain> E = mat_matrix(node_coords, r, s, t);
            StaticMatrix<n_strain, D * N> B = this->strain_displacements(node_coords, r, s, t, det);
            StaticMatrix<D * N, D * N> res = B.transpose() * (E * B) * det;
            return StaticMatrix<D * N, D * N>(res);
        };

    StaticMatrix<D * N, D * N> stiffness = integration_scheme().integrate(func);
    stiffness = 0.5 * (stiffness + stiffness.transpose()); // Symmetrize

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = stiffness;
    return mapped;
}

//-----------------------------------------------------------------------------
// mass
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::mass(NodeData& position, Precision* buffer) {
    logging::error(_material != nullptr, "no _material assigned to element ", elem_id);
    logging::error(_material->has_density(), "_material has no density assigned at element ", elem_id);

    Precision density = _material->get_density();

    StaticMatrix<N, D> node_coords = this->node_coords_global(position);

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, node_coords, density](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det;
            StaticMatrix<D, D> jac = this->jacobian(node_coords, r, s, t);
            StaticMatrix<N, N> shape_func_mass = this->shape_function(r, s, t) * this->shape_function(r, s, t).transpose();
            det = jac.determinant();

            // Expand the mass matrix from N x N to D * N x D * N
            StaticMatrix<D * N, D * N> mass_local = StaticMatrix<D * N, D * N>::Zero();

            for (Index i = 0; i < N; i++) {
                for (Index j = 0; j < N; j++) {
                    for (Dim d = 0; d < D; d++) {
                        mass_local(D * i + d, D * j + d) = shape_func_mass(i, j);
                    }
                }
            }

            return mass_local * det * density;
    };

    StaticMatrix<D * N, D * N> mass = integration_scheme().integrate(func);

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = mass;
    return mapped;
}


//-----------------------------------------------------------------------------
// dimensions
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::dimensions() {
    return D;
}

//-----------------------------------------------------------------------------
// n_nodes
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::n_nodes() {
    return node_ids.size();
}

//-----------------------------------------------------------------------------
// nodes
//-----------------------------------------------------------------------------
template<Index N>
ID*
SolidElement<N>::nodes() {
    return &node_ids[0];
}

//-----------------------------------------------------------------------------
// n_integration_points
//-----------------------------------------------------------------------------
template<Index N>
Dim
SolidElement<N>::n_integration_points() {
    return integration_scheme().count();
}

//-----------------------------------------------------------------------------
// volume
//-----------------------------------------------------------------------------
template<Index N>
Precision
SolidElement<N>::volume(NodeData& node_coords) {
    StaticMatrix<N, D> node_coords_glob = this->node_coords_global(node_coords);

    std::function<Precision(Precision, Precision, Precision)> func =
        [this, node_coords_glob](Precision r, Precision s, Precision t) -> Precision {
            Precision det = jacobian(node_coords_glob, r, s, t).determinant();
            return det;
        };

    Precision volume = integration_scheme().integrate(func);
    return volume;
}

//-----------------------------------------------------------------------------
// apply_vload
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::apply_vload(NodeData& node_coords, NodeData& node_loads, Vec3 load) {
    StaticMatrix<N, D> node_coords_glob = this->node_coords_global(node_coords);

    std::function<StaticMatrix<N, 1>(Precision, Precision, Precision)> func =
        [this, node_coords_glob](Precision r, Precision s, Precision t) -> StaticMatrix<N, 1> {
            Precision det = jacobian(node_coords_glob, r, s, t).determinant();
            StaticMatrix<N, 1> shape_func = shape_function(r, s, t);
            return (det * shape_func);
        };

    StaticMatrix<N, 1> nodal_impact = integration_scheme().integrate(func);

    ID local_id = 0;
    for (auto n_id : node_ids) {
        node_loads(n_id, 0) += nodal_impact(local_id) * load(0);
        node_loads(n_id, 1) += nodal_impact(local_id) * load(1);
        node_loads(n_id, 2) += nodal_impact(local_id) * load(2);
        local_id++;
    }
}

//-----------------------------------------------------------------------------
// apply tload
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::apply_tload(NodeData& node_coords, NodeData& node_loads, NodeData& node_temp, Precision ref_temp) {

    StaticMatrix<N, D> node_coords_glob = this->node_coords_global(node_coords);
    StaticMatrix<N, 1> node_temp_glob   = this->nodal_data<1>(node_temp);

    // adjust the temperature field to the reference temperature if its nan or inf
    for (Index i = 0; i < (Index)node_temp_glob.size(); i++) {
        if (std::isnan(node_temp_glob(i)) || std::isinf(node_temp_glob(i))) {
            node_temp_glob(i) = ref_temp;
        }
    }

    // error out if no _material is assigned
    logging::error(_material != nullptr, "no _material assigned to element ", elem_id);
    logging::error(_material->has_elasticity(), "_material has no elasticity components assigned at element ", elem_id);
    logging::error(_material->has_thermal_expansion(), "_material has no thermal expansion assigned at element ", elem_id);

    std::function<StaticMatrix<D*N,1>(Precision, Precision, Precision)> func =
        [this, node_coords_glob, &node_temp_glob, &ref_temp](Precision r, Precision s, Precision t) -> StaticMatrix<D*N,1> {

        // get temperature at integration point
        auto shape_func = shape_function(r, s, t);
        Precision temp = shape_func.dot(node_temp_glob);

        // create strain vector
        Precision strain_value = _material->get_thermal_expansion() * (temp - ref_temp);
        Vec6 strain{strain_value, strain_value, strain_value, 0, 0, 0};

        // stress tensor
        StaticMatrix<n_strain, n_strain> E = mat_matrix(node_coords_glob, r, s, t);
        auto stress = E * strain;

        // compute strain-displacement matrix
        Precision det;
        auto B = strain_displacements(node_coords_glob, r, s, t, det);

        auto res = B.transpose() * stress * det;
        return res;
    };

    StaticMatrix<D*N, 1> nodal_impact = integration_scheme().integrate(func);

    ID local_id = 0;
    for (auto n_id : node_ids) {
        node_loads(n_id, 0) += nodal_impact(local_id * 3 + 0);
        node_loads(n_id, 1) += nodal_impact(local_id * 3 + 1);
        node_loads(n_id, 2) += nodal_impact(local_id * 3 + 2);
        local_id++;
    }
}


//-----------------------------------------------------------------------------
// compute_stress_strain_nodal
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::compute_stress_strain_nodal(NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain) {
    auto local_node_coords = this->node_coords_local();
    auto global_node_coords = this->node_coords_global(node_coords);

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
            StaticMatrix<n_strain, n_strain> E = this->mat_matrix(global_node_coords, r, s, t);

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

            compute_stress_strain(node_coords, displacement, ip_stress, ip_strain, ip_xyz);
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
SolidElement<N>::compute_stress_strain(NodeData& node_coords, NodeData& displacement, NodeData& stress, NodeData& strain, NodeData& xyz) {
    auto global_node_coords = this->node_coords_global(node_coords);

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
        StaticMatrix<n_strain, n_strain> E = this->mat_matrix(global_node_coords, r, s, t);

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
SolidElement<N>::compute_compliance(NodeData& node_coords, NodeData& displacement, ElementData& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(node_coords, buffer);

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0) = strain_energy;
}

//-----------------------------------------------------------------------------
// test_implementation
//-----------------------------------------------------------------------------
template<Index N>
template<class ElementType>
bool
SolidElement<N>::test_implementation(bool print) {
    // Create an instance of the element type
    std::array<ID, N> nodeArray;
    for (size_t i = 0; i < N; i++) {
        nodeArray[i] = static_cast<ID>(i);    // Just initializing to example node IDs
    }

    ElementType el(0, nodeArray);

    // Test 1: Checking node values
    auto               node_coords = el.node_coords_local();
    StaticMatrix<N, N> globalMatrix;
    globalMatrix.setZero();

    for (size_t i = 0; i < N; i++) {
        Precision          r               = node_coords(i, 0);
        Precision          s               = node_coords(i, 1);
        Precision          t               = node_coords(i, 2);

        StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);
        for (size_t j = 0; j < N; j++) {
            globalMatrix(i, j) = shapeFuncValues(j);
        }
    }
    if (print)
        std::cout << globalMatrix << std::endl;

    // Test 2: Checking shape function sum
    const Precision step      = 0.2;
    const Precision tolerance = 1e-6;

    for (Precision r = -1; r <= 1; r += step) {
        for (Precision s = -1; s <= 1; s += step) {
            for (Precision t = -1; t <= 1; t += step) {
                StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);

                Precision          sum             = 0;
                for (size_t j = 0; j < N; j++) {
                    sum += shapeFuncValues(j);
                }

                if (std::abs(sum - 1.0) > tolerance) {
                    if (print)
                        std::cout << "Sum of shape functions at (r, s, t) = (" << r << ", " << s << ", " << t
                                  << ") is not 1. Actual sum: " << sum << std::endl;
                    return false;
                }
            }
        }
    }

    const Precision delta = 1e-6;
    for (Precision r = -1; r <= 1; r += step) {
        for (Precision s = -1; s <= 1; s += step) {
            for (Precision t = -1; t <= 1; t += step) {

                // Derivative from the function
                StaticMatrix<N, D> true_derivatives = el.shape_derivative(r, s, t);

                // Compute finite differences for each direction
                StaticMatrix<N, 1> shapeFuncValues_r_plus_delta = el.shape_function(r + delta, s, t);
                StaticMatrix<N, 1> shapeFuncValues_s_plus_delta = el.shape_function(r, s + delta, t);
                StaticMatrix<N, 1> shapeFuncValues_t_plus_delta = el.shape_function(r, s, t + delta);

                StaticMatrix<N, 1> shapeFuncValues_r_minu_delta = el.shape_function(r - delta, s, t);
                StaticMatrix<N, 1> shapeFuncValues_s_minu_delta = el.shape_function(r, s - delta, t);
                StaticMatrix<N, 1> shapeFuncValues_t_minu_delta = el.shape_function(r, s, t - delta);

                StaticMatrix<N, 1> shapeFuncValues = el.shape_function(r, s, t);

                StaticMatrix<N, D> finite_diff_derivatives;

                for (size_t j = 0; j < N; j++) {
                    finite_diff_derivatives(j, 0) = (shapeFuncValues_r_plus_delta(j) - shapeFuncValues_r_minu_delta(j)) / (2 * delta);  // dr
                    finite_diff_derivatives(j, 1) = (shapeFuncValues_s_plus_delta(j) - shapeFuncValues_s_minu_delta(j)) / (2 * delta);  // ds
                    finite_diff_derivatives(j, 2) = (shapeFuncValues_t_plus_delta(j) - shapeFuncValues_t_minu_delta(j)) / (2 * delta);  // dt
                }

                // Compare true derivatives with finite differences
                for (size_t j = 0; j < N; j++) {
                    for (size_t d = 0; d < D; d++) {
                        if (std::abs(true_derivatives(j, d) - finite_diff_derivatives(j, d)) > tolerance) {
                            if (print)
                                std::cout << "Mismatch in derivative at (r, s, t) = (" << r << ", " << s << ", " << t
                                          << ") in direction " << d
                                          << " from shape function " << j
                                          << ". True derivative: " << true_derivatives(j, d)
                                          << ", Finite Difference: " << finite_diff_derivatives(j, d) << std::endl;
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}

}  // namespace fem::model
