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
// strain_displacement
//-----------------------------------------------------------------------------
template<Index N>
StaticMatrix<SolidElement<N>::n_strain, SolidElement<N>::D * N>
    SolidElement<N>::strain_displacement(const StaticMatrix<N, D>& shape_der_global) {
    StaticMatrix<n_strain, D * N> B {};
    B.setZero();

    for (Index j = 0; j < N; j++) {
        Dim r1   = j * 3;
        Dim r2   = r1 + 1;
        Dim r3   = r1 + 2;

        B(0, r1) = shape_der_global(j, 0);    // dx/dr
        B(1, r2) = shape_der_global(j, 1);    // dy/ds
        B(2, r3) = shape_der_global(j, 2);    // dz/dt

        B(3, r2) = shape_der_global(j, 2);    // dy/dt
        B(3, r3) = shape_der_global(j, 1);    // dz/ds

        B(4, r1) = shape_der_global(j, 2);    // dz/dr
        B(4, r3) = shape_der_global(j, 0);    // dx/dt

        B(5, r1) = shape_der_global(j, 1);    // dx/ds
        B(5, r2) = shape_der_global(j, 0);    // dy/dt
    }

    return B;
}

template<Index N>
StaticMatrix<SolidElement<N>::n_strain, SolidElement<N>::n_strain>
    SolidElement<N>::material_matrix(Precision r, Precision s, Precision t) {

    logging::error(material() != nullptr, "no _material assigned to element ", elem_id);
    logging::error(material()->has_elasticity(), "_material has no elasticity components assigned at element ", elem_id);

    Precision scaling = 1;
    if (this->_model_data->elem_data.has(TOPO_STIFFNESS)) {
        scaling = this->_model_data->elem_data.get(TOPO_STIFFNESS)(this->elem_id);
    }

    if (this->_model_data->elem_data.has(TOPO_ANGLES)) {
        Vec3                   angles       = this->_model_data->elem_data.get(TOPO_ANGLES).row(this->elem_id);
        cos::RectangularSystem rot          = cos::RectangularSystem::euler(angles(0), angles(1), angles(2));

        Vec3                   point_global = this->interpolate<D>(this->node_coords_global(), r, s, t);
        Vec3                   point_local  = rot.to_local(point_global);

        auto result = this->material()->elasticity()->template get_transformed<D>(rot.get_axes(point_local));
        return scaling * result;
    } else {
        return scaling * this->material()->elasticity()->template get<D>();
    }
}

template<Index N>
template<Dim K>
StaticVector<K> SolidElement<N>::interpolate(StaticMatrix<N, K> data, Precision r, Precision s, Precision t) {
    StaticMatrix<N, 1> shape_func = shape_function(r, s, t);
    StaticVector<K> res {};
    for (Index i = 0; i < K; i++) {
        res(i) = shape_func.dot(data.col(i));
    }
    return res;
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
// stiffness
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::stiffness(Precision* buffer) {

    StaticMatrix<N, D> node_coords = this->node_coords_global();

    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &node_coords](Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N> {
            Precision det;
            StaticMatrix<n_strain, D * N> B = this->strain_displacements(node_coords, r, s, t, det);
            StaticMatrix<n_strain, n_strain> E = material_matrix(r,s,t);
            StaticMatrix<D * N, D * N> res = B.transpose() * (E * B) * det;
            return StaticMatrix<D * N, D * N>(res);
        };

    StaticMatrix<D * N, D * N> stiffness = integration_scheme().integrate(func);
    stiffness = 0.5 * (stiffness + stiffness.transpose()); // Symmetrize

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = stiffness;
    return mapped;
}

template<Index N>
MapMatrix
SolidElement<N>::stiffness_geom(Precision* buffer, IPData& ip_stress, int ip_start_idx) {
    StaticMatrix<N, D> node_coords = this->node_coords_global();

    Index ip_counter = 0;

    // Funktor OHNE w-Parameter; integrate() übernimmt die Gewichte.
    std::function<StaticMatrix<D * N, D * N>(Precision, Precision, Precision)> func =
        [this, &node_coords, &ip_stress, ip_start_idx, &ip_counter]
        (Precision r, Precision s, Precision t) -> StaticMatrix<D * N, D * N>
    {

        // ---- Stress at IP (Voigt = [xx, yy, zz, yz, zx, xy]) ----
        auto v = ip_stress.row(ip_start_idx + ip_counter++);
        Eigen::Matrix<Precision,3,3> sigma;
        sigma << v(0), v(5), v(4),
            v(5), v(1), v(3),
            v(4), v(3), v(2);

        // ---- Local -> global shape derivatives, Jacobian, det ----
        StaticMatrix<N, D> dN_local = this->shape_derivative(r, s, t);
        StaticMatrix<D, D> J        = this->jacobian(node_coords, r, s, t);
        StaticMatrix<D, D> Jinv     = J.inverse();
        Precision det               = J.determinant();
        StaticMatrix<N, D> dN       = (Jinv * dN_local.transpose()).transpose(); // global grads per node

        // ---- Geometric stiffness integrand (scalar s_ab * I3) ----
        StaticMatrix<D * N, D * N> Kg = StaticMatrix<D * N, D * N>::Zero();

        for (Index a = 0; a < N; ++a)
        {
            Eigen::Matrix<Precision,3,1> gNa;
            gNa << dN(a,0), dN(a,1), dN(a,2);

            for (Index b = 0; b < N; ++b)
            {
                Eigen::Matrix<Precision,3,1> gNb;
                gNb << dN(b,0), dN(b,1), dN(b,2);

                Precision s_ab = (gNa.transpose() * sigma * gNb)(0,0); // scalar

                // add s_ab * I3 to block (a,b)
                for (int d = 0; d < 3; ++d)
                    Kg(3*a + d, 3*b + d) += s_ab * det;
            }
        }

        return Kg;

    };

    StaticMatrix<D * N, D * N> Kg = integration_scheme().integrate(func);
    Kg = 0.5 * (Kg + Kg.transpose()); // Numerisch symmetrisieren

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = Kg;
    return mapped;
}

//-----------------------------------------------------------------------------
// mass
//-----------------------------------------------------------------------------
template<Index N>
MapMatrix
SolidElement<N>::mass(Precision* buffer) {
    logging::error(material() != nullptr, "no material assigned to element ", elem_id);
    logging::error(material()->has_density(), "material has no density assigned at element ", elem_id);

    Precision density = material()->get_density();

    StaticMatrix<N, D> node_coords = this->node_coords_global();

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
// volume
//-----------------------------------------------------------------------------
template<Index N>
Precision
SolidElement<N>::volume() {
    StaticMatrix<N, D> node_coords_glob = this->node_coords_global();

    std::function<Precision(Precision, Precision, Precision)> func =
        [this, node_coords_glob](Precision r, Precision s, Precision t) -> Precision {
            Precision det = jacobian(node_coords_glob, r, s, t).determinant();
            return det;
        };

    Precision volume = integration_scheme().integrate(func);
    return volume;
}



}  // namespace fem::model
