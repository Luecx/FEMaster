/**
 * @file element_solid_.ipp
 * @brief Implements solid topology metadata and natural recovery coordinates.
 *
 * These adapters expose translational DOFs, connectivity and constitutive
 * point counts through ElementInterface. Recovery locations use topology node
 * coordinates or the stiffness quadrature; geometry and material assembly
 * remain in the other SolidElement implementation files.
 *
 * @see SolidElement
 */

#pragma once

namespace fem::model {

template<Index N>
ElDofs SolidElement<N>::dofs() const {
    return ElDofs {true, true, true, false, false, false};
}

template<Index N>
Dim SolidElement<N>::dimensions() const {
    return D;
}

template<Index N>
Dim SolidElement<N>::n_nodes() const {
    return node_ids.size();
}

template<Index N>
const ID* SolidElement<N>::nodes() const {
    return &node_ids[0];
}

template<Index N>
Dim SolidElement<N>::num_ip() const {
    return integration_scheme_stiffness().count();
}

/**
 * Provides natural nodal coordinates for element-nodal result recovery.
 *
 * Rows follow the element connectivity and columns contain r, s and t. Concrete
 * reduced formulations may override these requested recovery locations.
 *
 * @return N-by-three matrix of natural output coordinates.
 */
template<Index N>
RowMatrix SolidElement<N>::stress_strain_nodal_rst() {
    // Preserve connectivity order when exposing natural nodal recovery locations
    auto local = this->node_coords_local();
    RowMatrix rst(static_cast<Index>(N), 3);
    for (Index i = 0; i < N; ++i) {
        rst(i, 0) = local(i, 0);
        rst(i, 1) = local(i, 1);
        rst(i, 2) = local(i, 2);
    }
    return rst;
}

/**
 * Provides natural coordinates of the constitutive integration points.
 *
 * The stiffness integration rule defines both these locations and the ordering
 * of material-state rows used by mechanical assembly and result recovery.
 *
 * @return One r, s, t row per constitutive integration point.
 */
template<Index N>
RowMatrix SolidElement<N>::stress_strain_ip_rst() {
    // Use constitutive quadrature order to match integration-point state storage
    const auto& scheme = this->integration_scheme_stiffness();
    RowMatrix rst(scheme.count(), 3);
    for (Index i = 0; i < scheme.count(); ++i) {
        rst(i, 0) = scheme.get_point(i).r;
        rst(i, 1) = scheme.get_point(i).s;
        rst(i, 2) = scheme.get_point(i).t;
    }
    return rst;
}

}  // namespace fem::model
