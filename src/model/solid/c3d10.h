/**
 * @file c3d10.h
 * @brief Declares the ten-node quadratic tetrahedral solid element.
 */

#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D10 : public SolidElement<10>{
    C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds);

    // Recreate only the persistent element topology. Instance-specific ids,
    // section assignment and runtime state are established by Model::compile().
    ElementPtr copy() const override { return std::make_shared<C3D10>(elem_id, node_ids); }

    std::string type_name() const override { return "C3D10"; }

    StaticMatrix<10, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<10, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<10, 3> node_coords_local() override;

    SurfacePtr surface(ID surface_id) override;

    const math::quadrature::Quadrature& integration_scheme() const override;
    const math::quadrature::Quadrature& integration_scheme_stiffness() const override;
};
} }
