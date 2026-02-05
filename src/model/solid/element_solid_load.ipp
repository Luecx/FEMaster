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
// apply_vload
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::apply_vload(Field& node_loads, Vec3 load) {
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
SolidElement<N>::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {

    StaticMatrix<N, D> node_coords_glob = this->node_coords_global();
    logging::error(node_temp.domain == FieldDomain::NODE,
                   "Temperature field '", node_temp.name, "' must be a node field");
    logging::error(node_temp.components == 1,
                   "Temperature field '", node_temp.name, "' must have 1 component");
    StaticMatrix<N, 1> node_temp_glob {};
    for (Index i = 0; i < N; ++i) {
        const Index row = static_cast<Index>(node_ids[i]);
        node_temp_glob(i) = node_temp(row, 0);
    }

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


//-----------------------------------------------------------------------------
// integrate_vec_field
//-----------------------------------------------------------------------------
template<Index N>
void
SolidElement<N>::integrate_vec_field(Field& node_loads,
                                     bool scale_by_density,
                                     const VecField& field)
{
    // Node coordinates
    StaticMatrix<N, D> node_coords_glob = this->node_coords_global();

    // Optional density scaling
    Precision rho = 1.0;
    if (scale_by_density) {
        auto mat = this->material();
        logging::error(mat != nullptr && mat->has_density(),
                       "SolidElement: material density is required when scale_by_density=true for element ", this->elem_id);
        rho = mat->get_density();
    }

    const auto& scheme = this->integration_scheme();
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto pt = scheme.get_point(ip);
        const Precision r = pt.r;
        const Precision s = pt.s;
        const Precision t = pt.t;
        const Precision w = pt.w;

        StaticMatrix<N, 1> Nvals = this->shape_function(r, s, t);
        const StaticMatrix<D, D> J = this->jacobian(node_coords_glob, r, s, t);
        const Precision detJ = J.determinant();

        Vec3 x_ip = Vec3::Zero();
        for (Index i = 0; i < N; ++i) x_ip += Nvals(i) * node_coords_glob.row(i);

        Vec3 f_ip = field(x_ip) * (rho * w * detJ);

        for (Index i = 0; i < N; ++i) {
            const ID n_id = this->node_ids[i];
            const Precision a = Nvals(i);
            node_loads(n_id, 0) += a * f_ip(0);
            node_loads(n_id, 1) += a * f_ip(1);
            node_loads(n_id, 2) += a * f_ip(2);
        }
    }
}


}  // namespace fem::model
