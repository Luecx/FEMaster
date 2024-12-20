#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D20R : public SolidElement<20>{
    C3D20R(ID p_elem_id, const std::array<ID, 20>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override;

    const quadrature::Quadrature& integration_scheme_mass() override;

    SurfacePtr                    surface(ID surface_id) override;

    StaticMatrix<20, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<20, 3> node_coords_local() override;

};

} }
