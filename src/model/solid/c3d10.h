#pragma once

#include "element_solid.h"

namespace fem { namespace model {

struct C3D10 : public SolidElement<10>{

    C3D10(ID pElemId, const std::array<ID, 10>& pNodeIds);
    std::string type_name() const override { return "C3D10"; }

    StaticMatrix<10, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<10, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<10, 3> node_coords_local() override;

    SurfacePtr                    surface(ID surface_id) override;

    const quadrature::Quadrature& integration_scheme() override;

    const quadrature::Quadrature& integration_scheme_mass() override;
};

} }
