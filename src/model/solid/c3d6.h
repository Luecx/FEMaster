/**
 * @file c3d6.h
 * @brief Declares the six-node linear wedge solid element.
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D6 : public SolidElement<6>{
    C3D6(ID p_elem_id, const std::array<ID, 6>& p_node_ids);

    // Recreate only the persistent element topology. Instance-specific ids and
    // runtime bindings are assigned by Model::compile() afterwards.
    ElementPtr copy() const override { return std::make_shared<C3D6>(elem_id, node_ids); }

    const math::quadrature::Quadrature& integration_scheme() const override;

    SurfacePtr surface(ID surface_id) override;

    StaticMatrix<6, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> node_coords_local() override;

    std::string type_name() const override { return "C3D6"; }
};
} }
