#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D13 : public SolidElement<13>{
    C3D13(ID p_elem_id, const std::array<ID, 13>& p_node_ids);
    std::string type_name() const override { return "C3D13"; }

    const quadrature::Quadrature& integration_scheme() override;

    SurfacePtr                    surface(ID surface_id) override;

    StaticMatrix<13, 1> shape_function(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> shape_derivative(Precision r, Precision s, Precision t) override;
    StaticMatrix<13, 3> node_coords_local() override;

};

} }
