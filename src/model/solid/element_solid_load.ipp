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

namespace fem::model {


//-----------------------------------------------------------------------------
// apply_vload
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::apply_vload(NodeData& node_loads, Vec3 load) {
    StaticMatrix<N, D> node_coords_glob = this->node_coords_global();

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
SolidElement<N>::apply_tload(NodeData& node_loads, NodeData& node_temp, Precision ref_temp) {

    StaticMatrix<N, D> node_coords_glob = this->node_coords_global();
    StaticMatrix<N, 1> node_temp_glob   = this->nodal_data<1>(node_temp);

    // adjust the temperature field to the reference temperature if its nan or inf
    for (Index i = 0; i < (Index)node_temp_glob.size(); i++) {
        if (std::isnan(node_temp_glob(i)) || std::isinf(node_temp_glob(i))) {
            node_temp_glob(i) = ref_temp;
        }
    }

    // error out if no _material is assigned
    logging::error(material() != nullptr, "no material assigned to element ", elem_id);
    logging::error(material()->has_elasticity(), "material has no elasticity components assigned at element ", elem_id);
    logging::error(material()->has_thermal_expansion(), "material has no thermal expansion assigned at element ", elem_id);

    std::function<StaticMatrix<D*N,1>(Precision, Precision, Precision)> func =
        [this, node_coords_glob, &node_temp_glob, &ref_temp](Precision r, Precision s, Precision t) -> StaticMatrix<D*N,1> {

        // get temperature at integration point
        auto shape_func = shape_function(r, s, t);
        Precision temp = shape_func.dot(node_temp_glob);

        // create strain vector
        Precision strain_value = material()->get_thermal_expansion() * (temp - ref_temp);
        Vec6 strain{strain_value, strain_value, strain_value, 0, 0, 0};

        // stress tensor
        auto mat_matrix = material_matrix(r,s,t);
        auto stress = mat_matrix * strain;

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


}  // namespace fem::model
