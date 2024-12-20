#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D6 : public SolidElement<6>{
    C3D6(ID p_elem_id, const std::array<ID, 6>& p_node_ids);

    const quadrature::Quadrature& integration_scheme() override;

    SurfacePtr                    surface(ID surface_id) override;

    StaticMatrix<6, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<6, 3> node_coords_local() override;

};

} }
