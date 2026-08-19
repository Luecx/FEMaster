/**
 * @file c3d20.h
 * @brief Declares the twenty-node quadratic hexahedral solid element.
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D20 : public SolidElement<20>{
    C3D20(ID p_elem_id, const std::array<ID, 20>& p_node_ids);

    // Recreate only the persistent element topology. Dense assembly ids and
    // runtime bindings are assigned after cloning by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D20>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D20"; }

    const math::quadrature::Quadrature& integration_scheme() const override;

    SurfacePtr surface(ID surface_id) override;

    StaticMatrix<20, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> node_coords_local() override;
};
} }
