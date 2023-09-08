#pragma once

#include "element.h"
#include "element_solid.h"

namespace fem { namespace model {

struct C3D4 : public SolidElement<4>{

    C3D4(ID pElemId, const std::array<ID, 4>& pNodeIds);

    StaticMatrix<4, 1> shape_function(Precision r, Precision s, Precision t) override;

    StaticMatrix<4, 3> shape_derivative(Precision r, Precision s, Precision t) override;

    StaticMatrix<4, 3> node_coords_local() override;

    const quadrature::Quadrature& integration_scheme() override {
        const static quadrature::Quadrature quad{quadrature::DOMAIN_ISO_TET, quadrature::ORDER_LINEAR};
        return quad;
    }
};

} }
