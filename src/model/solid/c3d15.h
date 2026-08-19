/**
 * @file c3d15.h
 * @brief Declares the fifteen-node quadratic wedge solid element.
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D15 : public SolidElement<15>{
    C3D15(ID p_elem_id, const std::array<ID, 15>& p_node_ids);

    // Recreate only the persistent element topology. Dense assembly ids and
    // runtime bindings are assigned after cloning by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D15>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D15"; }

    const math::quadrature::Quadrature& integration_scheme() const override;

    SurfacePtr surface(ID surface_id) override;

    StaticMatrix<15, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<15, 3> node_coords_local() override;
};
} }
