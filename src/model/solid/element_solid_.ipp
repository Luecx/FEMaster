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

namespace fem::model {


//-----------------------------------------------------------------------------
// dofs
//-----------------------------------------------------------------------------
template<Index N>
ElDofs
SolidElement<N>::dofs() {
    return ElDofs {true, true, true, false, false, false};
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

template<Index N>
RowMatrix SolidElement<N>::stress_strain_nodal_rst() {
    auto local = this->node_coords_local();
    RowMatrix rst(static_cast<Index>(N), 3);
    for (Index i = 0; i < N; ++i) {
        rst(i, 0) = local(i, 0);
        rst(i, 1) = local(i, 1);
        rst(i, 2) = local(i, 2);
    }
    return rst;
}

template<Index N>
RowMatrix SolidElement<N>::stress_strain_ip_rst() {
    const auto& scheme = this->integration_scheme();
    RowMatrix rst(scheme.count(), 3);
    for (Index i = 0; i < scheme.count(); ++i) {
        rst(i, 0) = scheme.get_point(i).r;
        rst(i, 1) = scheme.get_point(i).s;
        rst(i, 2) = scheme.get_point(i).t;
    }
    return rst;
}


}  // namespace fem::model
